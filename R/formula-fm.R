# Formula helpers for two-field factorization machines

#' Set up a Factorization-Machine Term in \pkg{brms}
#'
#' Set up a second-order factorization machine for two categorical fields.
#' The function does not evaluate its arguments; it exists purely to describe
#' a predictor term to the formula parser.
#'
#' @param field1,field2 Two untransformed variables identifying the levels of
#'   the fields whose interaction is to be factorized. Factors, character or
#'   logical vectors, and integer-valued identifiers are supported.
#' @param k Positive integer giving the effective rank of the interaction.
#' @param main Logical. If \code{TRUE} (the default), include a separate
#'   sum-to-zero main effect for each field. Set this to \code{FALSE} when the
#'   corresponding main effects are supplied elsewhere in the formula.
#'
#' @details
#' For observation \eqn{n}, with field indices \eqn{g_n} and \eqn{h_n}, the
#' term contributes
#' \deqn{a_{g_n} + b_{h_n} +
#' s\sqrt{G H}\,q_{g_n}^{T}\mathop{diag}(d)r_{h_n}}
#' to the linear predictor when \code{main = TRUE}; the first two terms are
#' omitted otherwise. The vectors \eqn{a} and \eqn{b} use Stan's native
#' \code{sum_to_zero_vector} type. Here \eqn{G} and \eqn{H} are the numbers
#' of levels, and \eqn{q} and \eqn{r} are centered semi-orthogonal frames:
#' their columns are orthonormal and each column sums to zero. They are formed
#' by Householder reflectors in the \eqn{G-1} and \eqn{H-1} dimensional
#' centered subspaces. The entries of \eqn{d} are strictly positive,
#' descending, and have squared sum one. Thus the interaction has zero row and
#' column sums and rank \code{k} almost surely.
#'
#' The reflector coordinates have fixed independent standard-normal priors,
#' which induce Haar-uniform frames. Each of the \code{k} reflectors per field
#' retains one prior-only radial auxiliary coordinate. A fixed uniform simplex
#' prior on ordered squared-spectrum gaps induces \eqn{d}, which is reported
#' in fitted models with the \code{sifm} parameter prefix. Class \code{sdfm} is
#' exactly the
#' root-mean-square interaction over the full \eqn{G} by \eqn{H} table, while
#' class \code{sdfm_main} controls the marginal scales of the two main effects.
#' Main-effect priors can be selected with the \code{group} argument of
#' \code{\link{set_prior}}. Ordering removes factor permutations and continuous
#' rotations almost surely, but paired column-sign reflections remain; the
#' interaction surface, rather than individual frame signs, is meaningful.
#' At full rank, at most one square frame is restricted to determinant one;
#' the other frame retains both determinant signs, so this removes a redundant
#' orientation component without excluding any interaction surface.
#' If both centered margins are square at rank \code{k}, the remaining relative
#' determinant sign is a genuine disconnected component of the full-rank
#' interaction, selected by the sign of the unrestricted frame's final scalar
#' reflector coordinate. The transform is discontinuous at zero in that
#' coordinate. This is most visible for a two-by-two rank-one interaction;
#' initialize and compare chains in both determinant components when that sign
#' is scientifically uncertain.
#'
#' The Stan reflector implementation is adapted from Seth Axen's
#' \href{https://github.com/sethaxen/stan_semiorthogonal_transforms}{
#' \code{stan_semiorthogonal_transforms}} (MIT license). The construction is
#' based on Stewart (\doi{10.1137/0717034}) and Nirwan and Bertschinger
#' (\url{https://proceedings.mlr.press/v97/nirwan19a.html}).
#'
#' The native sum-to-zero main effects and tuple-based reflector transform
#' require CmdStan 2.36 or newer. Accordingly, models containing \code{fm}
#' terms currently require
#' \code{backend = "cmdstanr"}. Prediction supports new combinations of levels
#' observed during fitting, but not previously unseen field levels.
#' If ordinary population-level main effects for either field are included,
#' set \code{main = FALSE}. Do not combine an \code{fm} term with the ordinary
#' full interaction of the same two fields.
#'
#' @return An object of class \code{fm_term}, which is interpreted by the
#'   formula parsing functions of \pkg{brms}.
#'
#' @examples
#' \dontrun{
#' # a rank-3 interaction plus both categorical main effects
#' fit <- brm(
#'   rating ~ fm(user, item, k = 3),
#'   data = ratings,
#'   backend = "cmdstanr"
#' )
#'
#' # use ordinary main effects and only factorize the interaction
#' fit_no_main <- brm(
#'   rating ~ user + item + fm(user, item, k = 3, main = FALSE),
#'   data = ratings,
#'   backend = "cmdstanr"
#' )
#' }
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{set_prior}}
#' @export
fm <- function(field1, field2, k = 5, main = TRUE) {
  fields <- as.list(substitute(list(field1, field2)))[-1L]
  fields <- ulapply(fields, deparse0, backtick = TRUE, width.cutoff = 500L)
  for (field in fields) {
    if (!nzchar(field) || !is_equal(all_vars(field), field)) {
      stop2("'fm' fields must be single untransformed variables.")
    }
    stopif_illegal_group(field)
  }
  if (identical(fields[1L], fields[2L])) {
    stop2("The two fields supplied to 'fm' must be different variables.")
  }
  k_numeric <- as_one_numeric(k)
  if (!is.finite(k_numeric) || !is_wholenumber(k_numeric) ||
      k_numeric < 1 || k_numeric > .Machine$integer.max) {
    stop2("Argument 'k' of 'fm' must be a positive integer.")
  }
  k <- as.integer(k_numeric)
  main <- as_one_logical(main)
  label <- deparse0(match.call())
  out <- list(
    term = fields, field1 = fields[1L], field2 = fields[2L],
    by = "NA", k = k, main = main, label = label
  )
  structure(out, class = "fm_term")
}

# Extract factorization-machine terms from a linear formula.
terms_fm <- function(formula) {
  out <- find_terms(formula, "fm")
  if (!length(out)) {
    return(NULL)
  }
  eterms <- lapply(out, eval2, envir = environment())
  fields <- unique(unlist(from_list(eterms, "term")))
  if (!length(fields)) {
    stop2("Two variables must be supplied to function 'fm'.")
  }
  out <- str2formula(out)
  attr(out, "allvars") <- str2formula(fields)
  out
}

# Gather data-dependent metadata for factorization-machine terms.
frame_fm <- function(x, data, basis = NULL) {
  if (is.formula(x)) {
    x <- brmsterms(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["fm"]]
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  fm_terms <- all_terms(form)
  if (length(basis) && length(basis) != length(fm_terms)) {
    stop2("Stored factorization-machine metadata does not match the formula.")
  }
  out <- data.frame(
    term = fm_terms, label = NA_character_, field1 = NA_character_,
    field2 = NA_character_, group1 = NA_character_, group2 = NA_character_,
    k = NA_integer_, main = NA, stringsAsFactors = FALSE
  )
  out$levels1 <- out$levels2 <- vector("list", nrow(out))
  for (i in seq_rows(out)) {
    fmi <- eval2(out$term[i])
    stopifnot(inherits(fmi, "fm_term"))
    fe_terms <- all_terms(x[["fe"]])
    fe_vars <- lapply(fe_terms, function(term) unique(all_vars(term)))
    duplicate_main <- fmi$term[vapply(fmi$term, function(field) {
      any(vapply(fe_vars, function(vars) {
        length(vars) == 1L && identical(vars, field)
      }, logical(1)))
    }, logical(1))]
    if (fmi$main && length(duplicate_main)) {
      stop2(
        "Field '", duplicate_main[1L], "' already has an ordinary ",
        "population-level main effect. Set 'main = FALSE' in term '",
        out$term[i], "' or remove the duplicate main effect."
      )
    }
    duplicate_interaction <- any(vapply(seq_along(fe_terms), function(j) {
      expr <- str2lang(fe_terms[j])
      is_interaction <- is.call(expr) && identical(expr[[1L]], as.name(":"))
      is_interaction && length(fe_vars[[j]]) == 2L &&
        setequal(fe_vars[[j]], fmi$term)
    }, logical(1)))
    if (duplicate_interaction) {
      stop2(
        "The ordinary interaction between '", fmi$field1, "' and '",
        fmi$field2, "' duplicates the factorization-machine interaction."
      )
    }
    values1 <- get_fm_values(fmi$field1, data)
    values2 <- get_fm_values(fmi$field2, data)
    levels1 <- extract_levels(values1)
    levels2 <- extract_levels(values2)
    if (length(basis)) {
      old <- basis[[i]]
      if (!identical(c(fmi$field1, fmi$field2), old$fields) ||
          !identical(fmi$k, old$k) || !identical(fmi$main, old$main)) {
        stop2("Stored factorization-machine metadata does not match term '",
              out$term[i], "'.")
      }
      levels1 <- old$levels1
      levels2 <- old$levels2
    }
    nlevels <- c(length(levels1), length(levels2))
    if (any(nlevels < 2L)) {
      stop2("Each field in an 'fm' term must have at least two levels. ",
            "Error occurred for term '", out$term[i], "'.")
    }
    max_rank <- min(nlevels - 1L)
    if (fmi$k > max_rank) {
      stop2("Rank 'k' of term '", out$term[i], "' must not exceed ",
            max_rank, ", the smaller field dimension after centering.")
    }
    out$label[i] <- paste0(rename(fmi$field1), ":", rename(fmi$field2))
    out$field1[i] <- fmi$field1
    out$field2[i] <- fmi$field2
    out$group1[i] <- rename(fmi$field1)
    out$group2[i] <- rename(fmi$field2)
    out$k[i] <- fmi$k
    out$main[i] <- fmi$main
    out$levels1[[i]] <- levels1
    out$levels2[[i]] <- levels2
  }
  pair_keys <- vapply(seq_rows(out), function(i) {
    paste(sort(c(out$field1[i], out$field2[i])), collapse = "\r")
  }, character(1))
  if (anyDuplicated(pair_keys)) {
    stop2("Duplicated factorization-machine field pairs are not allowed.")
  }
  main_fields <- unlist(Map(
    function(x, y, main) if (main) c(x, y),
    out$field1, out$field2, out$main
  ))
  if (anyDuplicated(main_fields)) {
    duplicate <- unique(main_fields[duplicated(main_fields)])[1L]
    stop2("Field '", duplicate, "' has more than one factorization-machine ",
          "main effect. Set 'main = FALSE' on the additional term.")
  }
  class(out) <- c("fmframe", "data.frame")
  out
}

get_fm_values <- function(field, data) {
  values <- eval2(field, data)
  invalid <- !is.atomic(values) || length(dim(values)) > 1L ||
    length(values) != nrow(data) || anyNA(values)
  if (!invalid && is.numeric(values)) {
    invalid <- any(!is.finite(values)) || !all(is_wholenumber(values))
  }
  if (invalid) {
    stop2("Factorization-machine field '", field, "' must be a complete ",
          "categorical vector or an integer-valued identifier.")
  }
  values
}

is.fmframe <- function(x) {
  inherits(x, "fmframe")
}

# Extract the field variables from FM formulas across model components.
get_fm_vars <- function(x) {
  collect <- function(formulas) {
    if (is.formula(formulas)) {
      return(all_vars(attr(formulas, "allvars") %||% formulas))
    }
    if (is.list(formulas)) {
      return(unique(unlist(lapply(formulas, collect), use.names = FALSE)))
    }
    character(0)
  }
  collect(get_effect(x, target = "fm"))
}

has_fm_terms <- function(x) {
  if (is.mvbrmsterms(x) || inherits(x, "mvbrmsframe")) {
    return(any(vapply(x$terms, has_fm_terms, logical(1))))
  }
  if (is.brmsterms(x) || inherits(x, "brmsframe")) {
    effects <- c(x$dpars, x$nlpars)
    return(any(vapply(effects, has_fm_terms, logical(1))))
  }
  is.btl(x) && is.formula(x[["fm"]])
}

fm_cmdstan_version_supported <- function(version) {
  if (length(version) != 1L || anyNA(version)) {
    return(FALSE)
  }
  comparison <- try(
    utils::compareVersion(as.character(version), "2.36.0"),
    silent = TRUE
  )
  !is_try_error(comparison) && comparison >= 0
}

validate_fm_backend <- function(x, backend) {
  if (!has_fm_terms(x)) {
    return(invisible(TRUE))
  }
  backend <- as_one_character(backend)
  if (backend == "mock") {
    return(invisible(TRUE))
  }
  if (backend != "cmdstanr") {
    stop2(
      "Factorization-machine terms use semi-orthogonal reflector frames ",
      "and currently require 'backend = \"cmdstanr\"' with CmdStan 2.36 ",
      "or newer."
    )
  }
  require_package("cmdstanr")
  version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  if (is_try_error(version) || !fm_cmdstan_version_supported(version)) {
    found <- if (is_try_error(version)) "not installed" else as.character(version)
    if (length(found) != 1L || anyNA(found)) {
      found <- "not installed"
    }
    stop2(
      "Factorization-machine terms require CmdStan 2.36 or newer; found ",
      found, "."
    )
  }
  invisible(TRUE)
}

# Preserve fitting levels for prediction with previously unseen pairs.
frame_basis_fm <- function(x, data, ...) {
  stopifnot(is.btl(x))
  fmframe <- if (is.bframel(x)) x$frame$fm else frame_fm(x, data)
  out <- vector("list", nrow(fmframe))
  for (i in seq_along(out)) {
    out[[i]] <- list(
      fields = c(fmframe$field1[i], fmframe$field2[i]),
      levels1 = fmframe$levels1[[i]], levels2 = fmframe$levels2[[i]],
      k = fmframe$k[i], main = fmframe$main[i]
    )
  }
  out
}

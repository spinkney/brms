# External tool output and user paths are literal data, not glue expressions.
.pnuts_stop <- function(...) {
  message <- paste0(..., collapse = "")
  stop2("{message}")
}

# Experimental native PNUTS backend. The sampler is installed separately.

.pnuts_tools <- function() {
  require_package("processx")
  require_package("jsonlite")
  resolve <- function(option, variable, command) {
    path <- getOption(option, Sys.getenv(variable, unset = command))
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
      .pnuts_stop("Set option '", option, "' or environment variable '", variable, "'.")
    }
    if (!file.exists(path)) path <- Sys.which(path)
    if (!nzchar(path) || !file.exists(path)) {
      .pnuts_stop("Cannot find ", command, ". Set option '", option, "' or '", variable, "'.")
    }
    # Preserve a virtual environment's Python symlink: resolving its final
    # target would silently select the base interpreter and lose its packages.
    file.path(normalizePath(dirname(path), mustWork = TRUE), basename(path))
  }
  list(
    executable = resolve("brms.pnuts.executable", "PNUTS_EXECUTABLE", "pnuts"),
    python = resolve("brms.pnuts.python", "PNUTS_PYTHON", "python3"),
    bridge = system.file("pnuts", "bridge.py", package = "brms", mustWork = TRUE),
    bridgestan_path = getOption("brms.pnuts.bridgestan", Sys.getenv("BRIDGESTAN"))
  )
}

.pnuts_call <- function(request, tools, directory, silent = 1, timeout = Inf) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  request_file <- tempfile("request-", tmpdir = directory, fileext = ".json")
  response_file <- tempfile("response-", tmpdir = directory, fileext = ".json")
  request$response_file <- response_file
  request$bridgestan_path <- tools$bridgestan_path
  jsonlite::write_json(request, request_file, auto_unbox = TRUE, digits = NA)
  result <- processx::run(
    tools$python, c(tools$bridge, request_file), error_on_status = FALSE,
    echo = silent == 0, timeout = timeout, cleanup_tree = TRUE
  )
  if (result$status != 0 || !file.exists(response_file)) {
    .pnuts_stop("PNUTS ", request$action, " failed.\n", result$stderr, result$stdout,
          "\nRequest and chain logs: ", directory)
  }
  jsonlite::read_json(response_file, simplifyVector = TRUE)
}

.parse_model_pnuts <- function(model, silent = 1, ...) {
  if (length(list(...))) .pnuts_stop("Unused arguments to the PNUTS parser.")
  tools <- .pnuts_tools()
  directory <- tempfile("brms-pnuts-parse-")
  dir.create(directory)
  file <- file.path(directory, "model.stan")
  writeLines(model, file)
  .pnuts_call(list(action = "parse", code = model, stan_file = file),
              tools, directory, silent = silent, timeout = 120)$code
}

.compile_model_pnuts <- function(model, threads, opencl, silent = 1,
                                 cache_dir = getOption("brms.pnuts.cache_dir",
                                   file.path(tools::R_user_dir("brms", "cache"), "pnuts")),
                                 make_args = character(), stanc_args = character(), ...) {
  if (use_opencl(opencl)) .pnuts_stop("The PNUTS backend does not support OpenCL.")
  if (length(list(...))) .pnuts_stop("Unused PNUTS compilation arguments: ", paste(names(list(...)), collapse = ", "))
  tools <- .pnuts_tools()
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- normalizePath(cache_dir, mustWork = TRUE)
  if (silent < 2) message("Compiling Stan program for PNUTS (BridgeStan cache)...")
  out <- .pnuts_call(list(action = "compile", code = model, cache_dir = cache_dir,
                         make_args = I(make_args), stanc_args = I(stanc_args)),
                     tools, cache_dir, silent = silent, timeout = 900)
  out$compile_args <- list(cache_dir = cache_dir, make_args = make_args, stanc_args = stanc_args)
  class(out) <- "brms_pnuts_model"
  out
}

.pnuts_integer <- function(x, name, lower = 0) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > .Machine$integer.max) {
    .pnuts_stop("PNUTS '", name, "' must be an integer >= ", lower, ".")
  }
  as.integer(x)
}

.pnuts_control <- function(control) {
  defaults <- list(adapt = "diag", geometry_frame = "pilot", method = "endpoint",
                   target_fidelity = 0.9, gamma = 0.5, step_size = 0.1,
                   max_depth = 10L, max_memory_mb = 1024,
                   chain_timeout = 1800)
  allowed <- c(names(defaults), "min_size", "threshold", "probability_floor", "power",
               "min_step_size", "max_step_size", "metric_regularization", "metric_blend",
               "geometry_pilot", "pilot_iterations", "pilot_geometry", "pilot_rank",
               "pilot_probes", "feature_limit", "state_dependent", "geometry_diagnostics")
  if (is.null(control)) control <- list()
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control))))) {
    .pnuts_stop("PNUTS 'control' must be a uniquely named list.")
  }
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    .pnuts_stop("Unknown PNUTS control: ", paste(unknown, collapse = ", "),
          ". Use target_fidelity (not adapt_delta), and max_depth (not max_treedepth).")
  }
  defaults[names(control)] <- control
  out <- defaults
  for (name in intersect(c("adapt", "geometry_frame", "method", "geometry_pilot", "pilot_geometry"), names(out))) {
    value <- out[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      .pnuts_stop("PNUTS control '", name, "' must be a single string.")
    }
  }
  if (!out$adapt %in% c("none", "step-size", "diag", "diag-fisher", "pilot-only")) {
    .pnuts_stop("Invalid PNUTS adaptation method.")
  }
  for (name in setdiff(names(out), c("adapt", "geometry_frame", "method", "geometry_pilot", "pilot_geometry",
                                     "state_dependent", "geometry_diagnostics"))) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      .pnuts_stop("PNUTS control '", name, "' must be a finite number.")
    }
  }
  if (out$chain_timeout <= 0) .pnuts_stop("PNUTS 'chain_timeout' must be positive.")
  for (name in intersect(c("state_dependent", "geometry_diagnostics"), names(out))) {
    out[[name]] <- as_one_logical(out[[name]])
  }
  out
}

.pnuts_initializations <- function(init, chains) {
  if (is.null(init) || is_equal(init, "random")) {
    return(rep(list(list(kind = "random", radius = 2)), chains))
  }
  if (is_equal(init, "0")) init <- 0
  if (is.numeric(init) && length(init) == 1L && is.finite(init) && init >= 0) {
    return(rep(list(list(kind = if (init == 0) "zero" else "random", radius = init)), chains))
  }
  if (is.character(init) && length(init) == 1L) init <- match.fun(init)
  if (is.function(init)) {
    init <- lapply(seq_len(chains), function(i) {
      if ("chain_id" %in% names(formals(init))) init(chain_id = i) else init()
    })
  }
  if (!is.list(init) || length(init) != chains || !all(vapply(init, is.list, logical(1)))) {
    .pnuts_stop("PNUTS init must be 0, a random radius, a function, or one named parameter list per chain.")
  }
  lapply(init, function(x) {
    if (!length(x)) return(list(kind = "random", radius = 2))
    if (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
      .pnuts_stop("Each PNUTS parameter initialization must have unique names.")
    }
    list(kind = "constrained", radius = 2, values = x)
  })
}

.fit_model_pnuts <- function(model, sdata, algorithm, iter, warmup, thin,
                             chains, cores, threads, opencl, init, exclude,
                             seed, control, silent, future, output_dir = NULL,
                             save_warmup = FALSE, refresh = 0, ...) {
  if (algorithm != "sampling") .pnuts_stop("PNUTS currently supports algorithm = 'sampling' only.")
  if (use_opencl(opencl)) .pnuts_stop("The PNUTS backend does not support OpenCL.")
  if (future) .pnuts_stop("PNUTS uses native chain processes; use 'cores' instead of 'future'.")
  if (length(list(...))) .pnuts_stop("Unused PNUTS sampling arguments: ", paste(names(list(...)), collapse = ", "))
  iter <- .pnuts_integer(iter, "iter", 1)
  warmup <- .pnuts_integer(warmup, "warmup")
  thin <- .pnuts_integer(thin, "thin", 1)
  chains <- .pnuts_integer(chains, "chains", 1)
  cores <- min(chains, .pnuts_integer(cores, "cores", 1))
  if (warmup >= iter) .pnuts_stop("PNUTS requires iter > warmup.")
  control <- .pnuts_control(control)
  if (control$adapt != "none" && warmup < 20) .pnuts_stop("PNUTS adaptation needs at least 20 warmup iterations.")
  if (isNA(seed)) seed <- sample.int(.Machine$integer.max - 1L, 1L)
  seed <- .pnuts_integer(seed, "seed")
  seeds <- as.integer((as.double(seed) + 104729 * (seq_len(chains) - 1L)) %% .Machine$integer.max)
  tools <- .pnuts_tools()
  if (!file.exists(model$library)) .pnuts_stop("PNUTS model library is missing; recompile the model.")
  if (is.null(output_dir)) output_dir <- tempfile("brms-pnuts-")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  if (any(file.exists(file.path(output_dir, c("data.json", "run.json", "chain-1.csv", "chain-1.csv.partial"))))) {
    .pnuts_stop("PNUTS output_dir already contains a run. Use a fresh directory: ", output_dir)
  }
  data_file <- file.path(output_dir, "data.json")
  # brms standata uses array dimensions to distinguish Stan vectors of size one.
  class(sdata) <- "list"
  jsonlite::write_json(sdata, data_file, auto_unbox = TRUE, digits = NA, factor = "integer")
  arguments <- character()
  for (name in setdiff(names(control), "chain_timeout")) {
    flag <- paste0("--", gsub("_", "-", name, fixed = TRUE))
    if (is.logical(control[[name]])) {
      if (control[[name]]) arguments <- c(arguments, flag)
    } else {
      arguments <- c(arguments, flag, as.character(control[[name]]))
    }
  }
  if (silent < 2) message("Start sampling with PNUTS (", chains, " chains)")
  request <- list(action = "sample", library = model$library, data_file = data_file,
                  executable = tools$executable, model_seed = seed,
                  seeds = I(seeds), initializations = .pnuts_initializations(init, chains),
                  draws = iter - warmup, warmup = warmup, cores = cores,
                  threads = as.integer(threads$threads %||% 1), output_dir = output_dir,
                  sampler_args = I(arguments), save_warmup = as_one_logical(save_warmup),
                  chain_timeout = control$chain_timeout)
  result <- .pnuts_call(request, tools, output_dir, silent = silent)
  .pnuts_read_fit(result$files, model, exclude, iter, warmup, thin, seed,
                  control, result, save_warmup)
}

.pnuts_read_fit <- function(files, model, exclude, iter, warmup, thin, seed,
                            control, run, save_warmup) {
  columns <- NULL
  samples <- diagnostics <- vector("list", length(files))
  warm <- warm_diag <- vector("list", length(files))
  for (i in seq_along(files)) {
    if (grepl("[.]partial$", files[i])) .pnuts_stop("Refusing partial PNUTS output.")
    text <- readLines(files[i], warn = FALSE)
    if (!any(grepl("^# seconds_total=", text)) || any(grepl("^# FAILED:", text))) {
      .pnuts_stop("PNUTS output has no successful completion record: ", files[i])
    }
    x <- utils::read.csv(files[i], comment.char = "#", check.names = FALSE)
    if (!is.null(columns) && !identical(columns, names(x))) .pnuts_stop("PNUTS chain columns differ.")
    columns <- names(x)
    expected <- iter - warmup + if (save_warmup) warmup else 0L
    if (nrow(x) != expected || !all(vapply(x, is.numeric, logical(1))) || !all(is.finite(as.matrix(x)))) {
      .pnuts_stop("Malformed or nonfinite PNUTS output: ", files[i])
    }
    if (!all(c("lp__", "warmup__", "stepsize__", "adapting__", "metric_updated__", "probe_grad_evals__") %in% names(x))) {
      .pnuts_stop("Missing required PNUTS output columns.")
    }
    keep <- which(x$warmup__ == 0)
    if (length(keep) != iter - warmup || length(unique(x$stepsize__[keep])) != 1L ||
        any(as.matrix(x[keep, c("adapting__", "metric_updated__", "probe_grad_evals__")]) != 0)) {
      .pnuts_stop("PNUTS retained draws do not have frozen tuning.")
    }
    diagnostic_names <- names(x)[grepl("__$", names(x))]
    draw_names <- names(x)[!grepl("__$", names(x)) | names(x) == "lp__"]
    draw_names <- draw_names[!sub("[.].*", "", draw_names) %in% exclude]
    draws <- x[draw_names]
    names(draws) <- repair_variable_names(names(draws))
    keep <- keep[seq.int(1L, length(keep), by = thin)]
    discarded <- which(x$warmup__ == 1)
    if (length(discarded)) discarded <- discarded[seq.int(1L, length(discarded), by = thin)]
    samples[[i]] <- as.matrix(draws[keep, , drop = FALSE])
    warm[[i]] <- as.matrix(draws[discarded, , drop = FALSE])
    diagnostics[[i]] <- as.matrix(x[keep, diagnostic_names, drop = FALSE])
    warm_diag[[i]] <- as.matrix(x[discarded, diagnostic_names, drop = FALSE])
  }
  as_array <- function(values) {
    out <- array(NA_real_, dim = c(nrow(values[[1]]), length(values), ncol(values[[1]])),
                 dimnames = list(NULL, NULL, colnames(values[[1]])))
    for (i in seq_along(values)) out[, i, ] <- values[[i]]
    posterior::as_draws_array(out)
  }
  draws <- as_array(samples)
  native <- as_array(diagnostics)
  flat_names <- colnames(samples[[1]])
  variables <- unique(sub("\\[.*", "", flat_names))
  time <- data.frame(warmup = rep(NA_real_, length(files)), sampling = NA_real_, total = run$statuses$seconds)
  meta <- list(model_name = "brms_pnuts", stan_variables = variables,
               variables = flat_names, stan_variable_sizes = .stanr_par_dims(flat_names, variables),
               method = "sampling", algorithm = "pnuts", engine = "pnuts", metric = control$geometry_frame,
               iter_warmup = warmup, iter_sampling = iter - warmup, thin = thin,
               save_warmup = save_warmup, num_chains = length(files), seed = seed,
               adapt_engaged = control$adapt != "none", max_treedepth = control$max_depth,
               target_fidelity = control$target_fidelity, gamma = control$gamma,
               step_size = vapply(diagnostics, function(x) x[1, "stepsize__"], numeric(1)),
               init = rep(NA_character_, length(files)), time = time,
               stan_version_major = NA_character_, stan_version_minor = NA_character_,
               stan_version_patch = NA_character_, stanc_version = NA_character_)
  csfit <- list(metadata = meta, inv_metric = NULL, step_size = as.list(meta$step_size),
                warmup_draws = as_array(warm), post_warmup_draws = draws,
                warmup_sampler_diagnostics = as_array(warm_diag), post_warmup_sampler_diagnostics = native)
  out <- .stanfit_from_csfit(csfit, files = files, model = model,
                            algorithm = "sampling", model_attr = "PNUTSModel")
  for (i in seq_along(out@sim$samples)) {
    attr(out@sim$samples[[i]], "args")$control <- control
    out@stan_args[[i]]$control <- control
  }
  attr(out, "pnuts") <- list(diagnostics = native, warmup_diagnostics = as_array(warm_diag),
                             control = control, run = run, files = files,
                             frozen_kernels = Filter(file.exists, paste0(files, ".adaptation")))
  repair_stanfit(out)
}

#' Extract PNUTS sampler diagnostics
#'
#' Native tree, numerical-fidelity and reflection diagnostics from a PNUTS fit.
#' These are not HMC acceptance probabilities or divergence indicators.
#' @param object A \code{brmsfit} fitted with \code{engine = "pnuts"}.
#' @param inc_warmup Include saved warmup diagnostics.
#' @return A \code{posterior::draws_array} with one variable per native diagnostic.
#' @export
pnuts_diagnostics <- function(object, inc_warmup = FALSE) {
  contains_draws(object)
  require_backend("pnuts", object)
  info <- attr(object$fit, "pnuts")
  if (as_one_logical(inc_warmup) && posterior::niterations(info$warmup_diagnostics) > 0) {
    return(posterior::bind_draws(info$warmup_diagnostics, info$diagnostics, along = "iteration"))
  }
  info$diagnostics
}

library(cmdstanr)

script_dir <- file.path(
  "analysis", "fm-s2z-diagnostics", "benchmark", "expanded"
)
if (!file.exists(file.path(script_dir, "scenarios.R"))) {
  stop("Run this script from the brms repository root.")
}
source(file.path(script_dir, "scenarios.R"))

root <- Sys.getenv("FM_EXPANDED_DIR")
if (!nzchar(root)) {
  stop("Set FM_EXPANDED_DIR to a fresh writable benchmark directory.")
}
dir.create(root, recursive = TRUE, showWarnings = FALSE)

env_integer <- function(name, default) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) {
    return(as.integer(default))
  }
  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

select_values <- function(name, choices) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) {
    return(choices)
  }
  selected <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(selected, choices)
  if (length(unknown)) {
    stop("Unknown values in ", name, ": ", paste(unknown, collapse = ", "))
  }
  selected
}

chains <- env_integer("FM_EXPANDED_CHAINS", 4L)
warmup <- env_integer("FM_EXPANDED_WARMUP", 500L)
sampling <- env_integer("FM_EXPANDED_SAMPLING", 500L)
parallel_chains <- min(chains, env_integer("FM_EXPANDED_CORES", chains))

scenarios <- make_expanded_fm_scenarios()
scenario_ids <- select_values("FM_EXPANDED_SCENARIOS", names(scenarios))
model_ids <- select_values(
  "FM_EXPANDED_MODELS", c("centered_raw", "semiorthogonal")
)
saveRDS(scenarios, file.path(root, "scenarios.rds"))

scenario_table <- do.call(rbind, lapply(scenarios, function(x) {
  data.frame(
    scenario = x$id, purpose = x$purpose, design = x$design,
    N1 = x$n1, N2 = x$n2, K = x$k,
    observations = length(x$centered_raw_data$y),
    heldout_cells = length(x$test_cells),
    special1 = x$semi_data$special1,
    special2 = x$semi_data$special2,
    truth_rms = x$truth$rms,
    prior_rms_mean = x$rms_prior$mean,
    prior_rms_sd = x$rms_prior$sd,
    prior_rms_q05 = x$rms_prior$quantiles[1L],
    prior_rms_q50 = x$rms_prior$quantiles[2L],
    prior_rms_q95 = x$rms_prior$quantiles[3L]
  )
}))
write.csv(scenario_table, file.path(root, "scenarios.csv"), row.names = FALSE)

model_files <- c(
  centered_raw = file.path(script_dir, "stan_fm_centered_raw.stan"),
  semiorthogonal = file.path(script_dir, "stan_fm_semiorthogonal.stan")
)
dir.create(file.path(root, "executables"), showWarnings = FALSE)
models <- lapply(model_ids, function(model_id) {
  file <- model_files[[model_id]]
  cmdstan_model(
    file,
    exe_file = file.path(root, "executables", model_id),
    include_paths = normalizePath(file.path("inst", "chunks")),
    quiet = FALSE
  )
})
names(models) <- model_ids

run_one <- function(scenario, model_id) {
  output_dir <- file.path(root, "runs", scenario$id, model_id)
  manifest_file <- file.path(output_dir, "manifest.rds")
  if (file.exists(manifest_file)) {
    message("Skipping completed run: ", scenario$id, " / ", model_id)
    return(invisible(readRDS(manifest_file)))
  }
  if (dir.exists(output_dir) && length(list.files(output_dir))) {
    stop("Refusing to mix outputs in incomplete directory: ", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data <- if (model_id == "centered_raw") {
    scenario$centered_raw_data
  } else {
    scenario$semi_data
  }
  message("Sampling ", scenario$id, " / ", model_id)
  fit <- models[[model_id]]$sample(
    data = data,
    seed = scenario$seeds$sampler,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = warmup,
    iter_sampling = sampling,
    adapt_delta = 0.85,
    max_treedepth = 20,
    refresh = 100,
    save_warmup = FALSE,
    sig_figs = 12,
    output_dir = output_dir,
    output_basename = paste0(scenario$id, "_", model_id)
  )
  manifest <- list(
    scenario = scenario$id, model = model_id,
    csv_files = basename(normalizePath(fit$output_files())),
    chains = chains, warmup = warmup, sampling = sampling,
    adapt_delta = 0.85, max_treedepth = 20,
    seed = scenario$seeds$sampler,
    code_file = normalizePath(model_files[[model_id]]),
    code_md5 = unname(tools::md5sum(model_files[[model_id]])),
    transform_md5 = unname(tools::md5sum(
      file.path("inst", "chunks", "fun_semiorthogonal_fm.stan")
    )),
    cmdstan_version = as.character(cmdstan_version()),
    time = fit$time()
  )
  saveRDS(manifest, manifest_file)
  invisible(manifest)
}

for (i in seq_along(scenario_ids)) {
  scenario <- scenarios[[scenario_ids[i]]]
  order <- model_ids
  if (i %% 2L == 0L) {
    order <- rev(order)
  }
  for (model_id in order) {
    run_one(scenario, model_id)
  }
}

message("Expanded FM sampling complete.")

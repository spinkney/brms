% Generated from documentation in R/backends-pnuts.R
\name{pnuts_diagnostics}
\alias{pnuts_diagnostics}
\title{Extract PNUTS sampler diagnostics}
\usage{
pnuts_diagnostics(object, inc_warmup = FALSE)
}
\arguments{
\item{object}{A \code{brmsfit} fitted with \code{engine = "pnuts"}.}
\item{inc_warmup}{Include saved warmup diagnostics.}
}
\value{
A \code{posterior::draws_array} with one variable per native diagnostic.
}
\description{
Native tree, numerical-fidelity and reflection diagnostics from a PNUTS fit.
These are not HMC acceptance probabilities or divergence indicators.
}
\details{
The experimental PNUTS backend requires a separately installed PNUTS executable
and a Python environment with BridgeStan 2.9.x and NumPy. Set options
\code{brms.pnuts.executable} and \code{brms.pnuts.python}, or environment
variables \code{PNUTS_EXECUTABLE} and \code{PNUTS_PYTHON}. The optional
\code{brms.pnuts.bridgestan} option specifies the BridgeStan source directory.

PNUTS compiles the generated Stan program through BridgeStan, runs native
chain processes, and returns constrained draws as a standard \code{brmsfit}.
Use \code{engine = "pnuts"} or \code{backend = "pnuts"} in \code{brm}.
Only \code{algorithm = "sampling"} is supported.

Controls include \code{target_fidelity} (default 0.9), \code{max_depth}
(10), \code{gamma} (0.5), \code{adapt} (\code{"diag"}), and
\code{geometry_frame} (\code{"pilot"}). Tuning freezes after warmup.
The default per-chain timeout is 1800 seconds, changed with
\code{control = list(chain_timeout = ...)}. \code{adapt_delta} and
\code{max_treedepth} are not PNUTS controls.

An optional \code{output_dir} passed to \code{brm} preserves CSV files,
initial positions, logs, run metadata and frozen adaptation states. It must
be fresh for each fit. Without it a temporary directory is used. Fits saved
with \code{saveRDS} retain draws and diagnostics, but do not bundle the
compiled library or raw output files. See the branch's \file{PNUTS.md} for
installation, supported operations and the full control list.
}
\seealso{\code{\link{brm}}, \code{\link{control_params}}}

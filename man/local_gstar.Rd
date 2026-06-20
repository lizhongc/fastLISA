\name{local_gstar}
\alias{local_gstar}
\title{Local Getis-Ord G*}
\usage{
local_gstar(
  x,
  listw,
  nsim = 999L,
  iseed = NULL,
  p.value = 0.05,
  n.cores = 1L,
  moments = FALSE
)
}
\arguments{
\item{x}{Numeric vector of length \code{n}.}

\item{listw}{A \code{listw} object from \pkg{spdep}.}

\item{nsim}{Integer; number of permutations. Default \code{999L}. Must be at least 1.}

\item{iseed}{Integer seed for RNG, or \code{NULL}.}

\item{p.value}{Numeric significance cutoff. Default \code{0.05}.}

\item{n.cores}{Integer; number of OpenMP threads. Default \code{1L};
set higher to use multiple cores. Ignored on platforms without OpenMP.}

\item{moments}{Logical; if \code{TRUE}, append permutation-distribution
moments \code{E.G*i}, \code{Var.G*i}, \code{Skew.G*i}, and
\code{Kurt.G*i}. Default \code{FALSE}.}
}
\value{
A numeric matrix of class \code{c("localG", "matrix", "array")}
  with columns \code{G*i}, \code{Z.G*i}, and \code{Pr(folded) Sim}.
  When \code{moments = TRUE}, the permutation-moment columns are appended.
  It has the following attributes:
  \describe{
    \item{cluster}{A significance-filtered factor with levels
      \code{Not significant}, \code{High-High}, \code{Low-Low},
      \code{Undefined} and \code{Isolated}.}
    \item{gstari}{Logical flag set to \code{TRUE} indicating local G*.}
    \item{call}{The matched call.}
  }
}
\description{
\code{local_gstar} computes the Getis-Ord local \eqn{G^*_i} statistic, the
self-inclusive companion of \code{\link{local_g}}: observation \eqn{i} is treated
as its own neighbour (weight 1). With \eqn{m_i} valid neighbours, row-standardised
neighbour weights \eqn{w_{ij}}, and global total \eqn{S = \sum_k x_k},
\deqn{G^*_i = \frac{\left(\frac{\sum_{j \in N_i} w_{ij} x_j}{\sum_j w_{ij}}\right) m_i + x_i}{(m_i + 1)\, S},}
i.e. the average value over the focal unit and its neighbours divided by the
global sum. Large \eqn{G^*_i} flags a hot spot and small \eqn{G^*_i} a cold spot,
with the focal unit included.
}
\details{
Inference uses a conditional permutation test (\code{nsim} reps), with the focal
\eqn{x_i} held fixed while neighbour values are permuted. The pseudo p-value is
folded (two-tailed),
\deqn{p_i = \frac{\min(g,\ \mathrm{nsim} - g) + 1}{\mathrm{nsim} + 1},}
where \eqn{g} counts permutations with
\eqn{G_i^{*\,\mathrm{perm}} \ge G_i^{*\,\mathrm{obs}}}. The standardised score
\eqn{Z.G^*_i = (G^*_i - E_{\mathrm{perm}}) / \sqrt{\mathrm{Var}_{\mathrm{perm}}}}
uses the permutation moments; \code{Skew.G*i}/\code{Kurt.G*i} (when
\code{moments = TRUE}) follow the \pkg{e1071} type-3 convention.

Observations with a missing \code{x} value are labelled \code{Undefined} and
observations with no neighbours are labelled \code{Isolated}; both receive
\code{NA} for the p-value, Z-score and moments. The C backend re-seeds its
random number generator per observation, so results are identical for any
\code{n.cores}; \code{n.cores} is ignored when the package is built without
OpenMP.
}
\references{
Getis, A. and Ord, J. K. (1992) The Analysis of Spatial Association by Use of
Distance Statistics. \emph{Geographical Analysis} \bold{24}(3), 189--206.
\doi{10.1111/j.1538-4632.1992.tb00261.x}

Ord, J. K. and Getis, A. (1995) Local Spatial Autocorrelation Statistics:
Distributional Issues and an Application. \emph{Geographical Analysis}
\bold{27}(4), 286--306. \doi{10.1111/j.1538-4632.1995.tb00912.x}
}
\examples{
\dontshow{if (requireNamespace("spdep", quietly = TRUE)) withAutoprint(\{ # examplesIf}
lw <- spdep::nb2listw(spdep::cell2nb(7, 7))
x  <- as.numeric(seq_len(49))
res <- local_gstar(x, lw, nsim = 99L, n.cores = 1L)
head(res)
\dontshow{\}) # examplesIf}
}

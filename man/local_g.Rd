\name{local_g}
\alias{local_g}
\title{Local Getis-Ord G}
\usage{
local_g(
  x,
  listw,
  nsim = 999L,
  iseed = NULL,
  p.value = 0.05,
  n.cores = 1L,
  moments = FALSE,
  p.method = c("count", "rank")
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
moments \code{E.Gi}, \code{Var.Gi}, \code{Skew.Gi}, and \code{Kurt.Gi}.
Default \code{FALSE}.}

\item{p.method}{Character; how the observed statistic is located within its
permutation distribution. \code{"count"} (default) is the standard rule: it
counts the permutations at least as extreme as the observed value, matching the
folded pseudo p-value reported by \pkg{rgeoda} and \pkg{spdep}. \code{"rank"} is
the ties-averaged alternative, using \pkg{spdep}'s averaged rank of the observed
value. The two differ only when a permuted value exactly ties the observed one
(discrete or tie-prone data); both return a folded (smaller-tail) value.}
}
\value{
A numeric matrix of class \code{c("localG", "matrix", "array")}
  with columns \code{Gi}, \code{Z.Gi}, and \code{Pr(folded) Sim}.
  When \code{moments = TRUE}, the permutation-moment columns are appended.
  It has the following attributes:
  \describe{
    \item{cluster}{A significance-filtered factor with levels
      \code{Not significant}, \code{High-High}, \code{Low-Low},
      \code{Undefined} and \code{Isolated}.}
    \item{gstari}{Logical flag set to \code{FALSE} indicating local G (not G*).}
    \item{call}{The matched call.}
  }
}
\description{
\code{local_g} computes the Getis-Ord local \eqn{G_i} statistic, a Local
Indicator of Spatial Association that detects local clustering of high values
(\dQuote{hot spots}) and low values (\dQuote{cold spots}). For observation
\eqn{i}, with row-standardised spatial weights \eqn{w^*_{ij}} and the focal
value \emph{excluded} from both the lag and the denominator,
\deqn{G_i = \frac{\sum_{j \ne i} w^*_{ij} x_j}{\sum_k x_k - x_i}.}
A large \eqn{G_i} indicates that \eqn{i} is surrounded by high values; a small
\eqn{G_i} indicates a low-value neighbourhood. \eqn{G_i} contains no self term;
see \code{\link{local_gstar}} for the self-inclusive \eqn{G^*_i}.
}
\details{
Inference uses a conditional permutation test: the focal value \eqn{x_i} is held
fixed while the neighbouring values are randomly permuted \code{nsim} times. The
pseudo p-value is folded (two-tailed),
\deqn{p_i = \frac{\min(g,\ \mathrm{nsim} - g) + 1}{\mathrm{nsim} + 1},}
where \eqn{g} is the number of permutations with
\eqn{G_i^{\mathrm{perm}} \ge G_i^{\mathrm{obs}}}. The standardised score is
\eqn{Z.G_i = (G_i - E_{\mathrm{perm}}) / \sqrt{\mathrm{Var}_{\mathrm{perm}}}},
computed from the permutation mean and variance; \code{Skew.Gi} and
\code{Kurt.Gi} (when \code{moments = TRUE}) follow the \pkg{e1071} type-3
convention.

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
res <- local_g(x, lw, nsim = 99L, n.cores = 1L)
head(res)
\dontshow{\}) # examplesIf}
}

\name{local_multigeary}
\alias{local_multigeary}
\title{Multivariate Local Geary's C}
\usage{
local_multigeary(
  df,
  listw,
  nsim = 999L,
  scale = TRUE,
  iseed = NULL,
  p.value = 0.05,
  n.cores = 1L,
  moments = FALSE,
  p.method = c("count", "rank")
)
}
\arguments{
\item{df}{A data.frame or matrix with selected variables.}

\item{listw}{A \code{listw} object from \pkg{spdep}.}

\item{nsim}{Integer; number of permutations. Default \code{999L}. Must be at least 1.}

\item{scale}{Logical; if \code{TRUE} (default), standardise each variable in R.}

\item{iseed}{Integer seed for RNG, or \code{NULL}.}

\item{p.value}{Numeric significance cutoff. Default \code{0.05}.}

\item{n.cores}{Integer; number of OpenMP threads. Default \code{1L};
set higher to use multiple cores. Ignored on platforms without OpenMP.}

\item{moments}{Logical; if \code{TRUE}, append permutation-distribution
moments \code{E.Ci}, \code{Var.Ci}, \code{Skew.Ci}, and \code{Kurt.Ci}.
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
A numeric matrix of class \code{c("localC", "matrix", "array")}
  with columns \code{Ci}, \code{Z.Ci}, and \code{Pr Sim}. When
  \code{moments = TRUE}, the permutation-moment columns are appended.
  It has the following attributes:
  \describe{
    \item{cluster}{A significance-filtered factor with levels Not significant,
      Positive, Negative, Undefined, and Isolated.}
    \item{call}{The matched call.}
  }
}
\description{
\code{local_multigeary} computes the multivariate local Geary's \eqn{C_i}
(Anselin 2019), the average across \eqn{K} variables of the univariate
squared-difference statistic. On the sample (\eqn{n-1}) standardised variables
\eqn{z^1,\dots,z^K} (when \code{scale = TRUE}) with row-standardised weights
\eqn{w^*_{ij}},
\deqn{C_i = \frac{1}{K} \sum_{v=1}^{K} \sum_j w^*_{ij} (z^v_i - z^v_j)^2.}
A small \eqn{C_i} indicates that \eqn{i} is similar to its neighbours across all
variables (positive association); a large \eqn{C_i} indicates multivariate
dissimilarity (a spatial outlier).
}
\details{
Inference uses a one-tailed conditional permutation test (\code{nsim} reps); each
replicate applies the same permuted neighbour configuration to every variable.
The observed \eqn{C_i} is compared with the permutation mean to choose the tail,
and \eqn{p_i = (t + 1)/(\mathrm{nsim} + 1)}, where \eqn{t} counts permuted
statistics in that tail. The standardised score is
\eqn{Z.C_i = (C_i - E_{\mathrm{perm}})/\sqrt{\mathrm{Var}_{\mathrm{perm}}}};
\code{Skew.Ci}/\code{Kurt.Ci} (when \code{moments = TRUE}) follow the \pkg{e1071}
type-3 convention. Significant units are labelled \code{Positive} (similar) or
\code{Negative} (dissimilar).

Rows with any missing value are labelled \code{Undefined} and observations with
no neighbours are labelled \code{Isolated}; both receive \code{NA} for the
p-value, Z-score and moments. The C backend re-seeds its random number generator
per observation, so results are identical for any \code{n.cores}; \code{n.cores}
is ignored when the package is built without OpenMP.
}
\references{
Anselin, L. (2019) A Local Indicator of Multivariate Spatial Association:
Extending Geary's c. \emph{Geographical Analysis} \bold{51}(2), 133--150.
\doi{10.1111/gean.12164}
}
\examples{
\dontshow{if (requireNamespace("spdep", quietly = TRUE)) withAutoprint(\{ # examplesIf}
lw <- spdep::nb2listw(spdep::cell2nb(7, 7))
x  <- as.numeric(seq_len(49))
df <- cbind(x, rev(x))
res <- local_multigeary(df, lw, nsim = 99L, n.cores = 1L)
head(res)
\dontshow{\}) # examplesIf}
}

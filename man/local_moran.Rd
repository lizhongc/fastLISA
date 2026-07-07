\name{local_moran}
\alias{local_moran}
\title{Univariate Local Moran's I}
\usage{
local_moran(
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

\item{moments}{Logical; if \code{TRUE}, append the permutation-distribution
moments \code{E.Ii}, \code{Var.Ii}, \code{Skew.Ii}, and \code{Kurt.Ii}
to the result matrix. Default \code{FALSE}.}

\item{p.method}{Character; how the observed statistic is located within its
permutation distribution. \code{"count"} (default) is the standard rule: it
counts the permutations at least as extreme as the observed value, matching the
folded pseudo p-value reported by \pkg{rgeoda} and \pkg{spdep}. \code{"rank"} is
the ties-averaged alternative, using \pkg{spdep}'s averaged rank of the observed
value. The two differ only when a permuted value exactly ties the observed one
(discrete or tie-prone data); both return a folded (smaller-tail) value.}
}
\value{
A numeric matrix of class \code{c("localmoran", "matrix", "array")}
  with \code{n} rows and 3 columns by default:
  \describe{
    \item{Ii}{Observed univariate Moran statistic.}
    \item{Z.Ii}{Standardised Z-score computed from permutation moments.}
    \item{Pr(folded) Sim}{rgeoda-style folded empirical permutation p-value.}
  }
  When \code{moments = TRUE}, the permutation-distribution columns
  \code{E.Ii}, \code{Var.Ii}, \code{Skew.Ii}, and \code{Kurt.Ii} are
  appended.
  The matrix has the following attributes:
  \describe{
    \item{quadr}{Moran scatter-plot quadrant classification.}
    \item{cluster}{A significance-filtered factor with levels
      \code{Not significant}, \code{High-High}, \code{Low-Low}, \code{Low-High},
      \code{High-Low}, \code{Undefined} and \code{Isolated}.}
  }
}
\description{
\code{local_moran} computes the univariate local Moran's \eqn{I_i}, the classic
Local Indicator of Spatial Association measuring whether a value coincides with
its neighbours' average. On the sample (\eqn{n-1}) standardised variable \eqn{z}
with row-standardised weights \eqn{w^*_{ij}},
\deqn{I_i = z_i \sum_j w^*_{ij} z_j = z_i\, \mathrm{lag}(z)_i.}
A positive \eqn{I_i} indicates similarity to neighbours (a High-High or Low-Low
cluster); a negative \eqn{I_i} indicates a spatial outlier (High-Low or
Low-High). Standardisation uses the sample standard deviation (\eqn{n-1}
denominator).
}
\details{
This is the \eqn{x = y} special case of bivariate local Moran's I, so the
computation delegates to \code{\link{local_moran_bv}} (with \code{scale = TRUE})
and relabels the columns. Inference uses a conditional permutation test
(\code{nsim} reps) with the folded two-tailed pseudo p-value
\eqn{p_i = (\min(g, \mathrm{nsim}-g)+1)/(\mathrm{nsim}+1)}, where \eqn{g} is the
number of permutations with \eqn{I^{\mathrm{perm}} \ge I^{\mathrm{obs}}}. The
standardised score is
\eqn{Z.I_i = (I_i - E_{\mathrm{perm}})/\sqrt{\mathrm{Var}_{\mathrm{perm}}}} and
\code{Skew.Ii}/\code{Kurt.Ii} (when \code{moments = TRUE}) use the \pkg{e1071}
type-3 convention. NA observations are \code{Undefined} and neighbourless
observations are \code{Isolated} (both \code{NA} in p-value, Z-score and moments).
The C backend re-seeds its random number generator per observation, so results
are identical for any \code{n.cores}.
}
\references{
Anselin, L. (1995) Local Indicators of Spatial Association---LISA.
\emph{Geographical Analysis} \bold{27}(2), 93--115.
\doi{10.1111/j.1538-4632.1995.tb00338.x}
}
\examples{
\dontshow{if (requireNamespace("spdep", quietly = TRUE)) withAutoprint(\{ # examplesIf}
lw <- spdep::nb2listw(spdep::cell2nb(7, 7))
x  <- as.numeric(seq_len(49))
res <- local_moran(x, lw, nsim = 99L, n.cores = 1L, moments = TRUE)
head(res)
\dontshow{\}) # examplesIf}
}

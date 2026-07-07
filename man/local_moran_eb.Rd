\name{local_moran_eb}
\alias{local_moran_eb}
\title{Local Moran's I with Empirical Bayes (EB) Rate}
\usage{
local_moran_eb(
  event,
  base,
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
\item{event}{Numeric vector of events (e.g. case counts).}

\item{base}{Numeric vector of populations at risk.}

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
A numeric matrix of class \code{c("local_moran_eb", "matrix")} with
  \code{n} rows and 3 columns by default:
  \describe{
    \item{Ii}{Observed Local Moran statistic computed on EB-standardised rates.}
    \item{Z.Ii}{Standardised Z-score computed from permutation moments.}
    \item{Pr(folded) Sim}{Folded empirical permutation p-value.}
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
    \item{call}{The matched call.}
    \item{nsim}{Number of simulations used.}
  }
}
\description{
\code{local_moran_eb} computes local Moran's I on Empirical Bayes (EB)
variance-stabilised rates, for event-count data observed over a population at
risk. Raw rates \eqn{p_i = \mathrm{event}_i / \mathrm{base}_i} from small
populations are noisy; EB standardisation shrinks them toward the global rate
\eqn{b = \sum \mathrm{event} / \sum \mathrm{base}} using a variance component
\eqn{\hat{a}},
\deqn{z_i = \frac{p_i - b}{\sqrt{\hat{a} + b/\mathrm{base}_i}},}
and local Moran's I is then computed on \eqn{z}. This follows the GeoDa/libgeoda
EBLocalMoran formulation.
}
\details{
Two standardisations are applied and are not redundant: the EB rate
standardisation above stabilises the rate variance, and the usual sample
(\eqn{n-1}) z-score standardisation that local Moran's I requires is then applied
to the EB rates internally. Because univariate Moran's I on a standardised
variable equals the bivariate statistic with \eqn{x = y}, the permutation engine
is reused via \code{\link{local_moran_bv}}. Inference is a conditional
permutation test (\code{nsim} reps); the folded two-tailed pseudo p-value and the
score \eqn{Z.I_i = (I_i - E_{\mathrm{perm}})/\sqrt{\mathrm{Var}_{\mathrm{perm}}}}
are as in \code{\link{local_moran}}, and \code{Skew.Ii}/\code{Kurt.Ii} (when
\code{moments = TRUE}) use the \pkg{e1071} type-3 convention.

Observations with NA \code{event}/\code{base} or \code{base} \eqn{\le 0} are
labelled \code{Undefined}; observations with no neighbours are \code{Isolated};
both receive \code{NA} for the p-value, Z-score and moments. The C backend
re-seeds its random number generator per observation, so results are identical
for any \code{n.cores}.
}
\references{
\enc{Assunção}{Assuncao}, R. M. and Reis, E. A. (1999) A new proposal to adjust
Moran's I for population density. \emph{Statistics in Medicine} \bold{18}(16),
2147--2162.
\doi{10.1002/(SICI)1097-0258(19990830)18:16<2147::AID-SIM179>3.0.CO;2-I}

Anselin, L. (1995) Local Indicators of Spatial Association---LISA.
\emph{Geographical Analysis} \bold{27}(2), 93--115.
\doi{10.1111/j.1538-4632.1995.tb00338.x}
}
\examples{
\dontshow{if (requireNamespace("spdep", quietly = TRUE)) withAutoprint(\{ # examplesIf}
lw    <- spdep::nb2listw(spdep::cell2nb(7, 7))
event <- as.numeric(seq_len(49))
base  <- rep(100, 49)
res <- local_moran_eb(event, base, lw, nsim = 99L, n.cores = 1L)
head(res)
\dontshow{\}) # examplesIf}
}

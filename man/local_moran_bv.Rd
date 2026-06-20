\name{local_moran_bv}
\alias{local_moran_bv}
\title{Bivariate Local Moran's I}
\usage{
local_moran_bv(
  x,
  y,
  listw,
  nsim = 999L,
  scale = TRUE,
  iseed = NULL,
  p.value = 0.05,
  n.cores = 1L,
  moments = FALSE
)
}
\arguments{
\item{x}{Numeric vector of length \code{n}; the first variable.}

\item{y}{Numeric vector of length \code{n}; the second variable (lagged).}

\item{listw}{A \code{listw} object from \pkg{spdep} (any style: "W", "B",
"C", etc., including custom distance-decay weights). Observations with
no neighbours receive cluster code 6 (Isolated).}

\item{nsim}{Integer; number of permutations for the pseudo p-value.
Default \code{999L}.}

\item{scale}{Logical; if \code{TRUE} (default), \code{x} and \code{y}
are standardised in R (sample std dev). Set to
\code{FALSE} only if you have pre-standardised the data.}

\item{iseed}{Integer seed for the RNG, or \code{NULL} (default) to use
the package default (123456789). Passed as the \code{seed} argument to
the C backend.}

\item{p.value}{Numeric; observations with
\eqn{p > \text{p.value}} are recoded to cluster 0 (Not
significant). Default \code{0.05}.}

\item{n.cores}{Integer; number of OpenMP threads. Default \code{1L};
set higher to use multiple cores. Ignored on platforms without OpenMP.}

\item{moments}{Logical; if \code{TRUE}, append the permutation-distribution
moments \code{E.Ibvi}, \code{Var.Ibvi}, \code{Skew.Ibvi}, and
\code{Kurt.Ibvi} to the result matrix. Default \code{FALSE}.}
}
\value{
A numeric matrix of class \code{c("localmoran", "matrix", "array")}
  with \code{n} rows and 3 columns by default:
  \describe{
    \item{Ibvi}{Observed bivariate Moran statistic.}
    \item{Z.Ibvi}{Standardised Z-score computed from permutation moments.}
    \item{Pr(folded) Sim}{rgeoda-style folded empirical permutation p-value.}
  }
  When \code{moments = TRUE}, the permutation-distribution columns
  \code{E.Ibvi}, \code{Var.Ibvi}, \code{Skew.Ibvi}, and \code{Kurt.Ibvi}
  are appended.
  The matrix has the following attributes:
  \describe{
    \item{quadr}{A \code{data.frame} with three factor columns
      (\code{mean}, \code{median}, \code{pysal}) giving the Moran
      scatter-plot quadrant for each observation, computed on the original
      (unscaled) data scale.}
    \item{cluster}{A \code{factor} representing cluster classification
      (Not significant, High-High, Low-Low, Low-High, High-Low, Undefined, Isolated).}
    \item{call}{The matched call.}
  }
}
\description{
\code{local_moran_bv} computes the bivariate local Moran's \eqn{I_{bv,i}}, which
correlates a variable \eqn{x} at \eqn{i} with the spatial lag of a second
variable \eqn{y} over \eqn{i}'s neighbours. On the sample (\eqn{n-1}) standardised
variables \eqn{z_x, z_y} (when \code{scale = TRUE}) with row-standardised weights
\eqn{w^*_{ij}},
\deqn{I_{bv,i} = z_{x,i} \sum_j w^*_{ij} z_{y,j} = z_{x,i}\, \mathrm{lag}(z_y)_i.}
A positive \eqn{I_{bv,i}} means \eqn{i}'s \eqn{x} value coincides with high lagged
\eqn{y} nearby; a negative value indicates spatial mismatch. The univariate local
Moran's I is the special case \eqn{x = y} (see \code{\link{local_moran}}). The
backend is plain C with optional OpenMP parallelism.
}
\details{
Inference uses a conditional permutation test (\code{nsim} reps): the neighbour
\eqn{y}-values are permuted with the focal observation held fixed. See the
\emph{P-values} and \emph{Cluster codes} sections below for the folded p-value and
the cluster coding. The standardised score is
\eqn{Z.I_{bv,i} = (I_{bv,i} - E_{\mathrm{perm}})/\sqrt{\mathrm{Var}_{\mathrm{perm}}}};
\code{Skew.Ibvi}/\code{Kurt.Ibvi} (when \code{moments = TRUE}) follow the
\pkg{e1071} type-3 convention. NA observations are \code{Undefined} and
neighbourless observations are \code{Isolated} (both \code{NA} in p-value, Z-score
and moments). The C backend re-seeds its random number generator per observation,
so results are identical for any \code{n.cores}.
}
\section{Standardisation}{

Both \code{x} and \code{y} are standardised using the \strong{sample}
standard deviation (\eqn{n-1} denominator) in R before computing the
statistic, consistent with the blisa backend.
}

\section{P-values}{

Permutation p-values use the folded two-tailed formula matching rgeoda:
\deqn{p = (\min(\#\{perm \ge obs\},\, \#\{perm < obs\}) + 1) \,/\, (nsim + 1)}
No normal approximation is computed.
}

\section{Cluster codes}{

The returned \code{cluster} factor attribute is based on integer codes 0–6:
\tabular{ll}{
  0 \tab Not significant \cr
  1 \tab High-High \cr
  2 \tab Low-Low \cr
  3 \tab Low-High \cr
  4 \tab High-Low \cr
  5 \tab Undefined (NA input) \cr
  6 \tab Isolated (no neighbours) \cr
}
Codes 5 and 6 are preserved regardless of the significance cutoff.
}

\references{
Anselin, L. (1995) Local Indicators of Spatial Association---LISA.
\emph{Geographical Analysis} \bold{27}(2), 93--115.
\doi{10.1111/j.1538-4632.1995.tb00338.x}
}

\examples{
\dontshow{if (requireNamespace("spdep", quietly = TRUE)) withAutoprint(\{ # examplesIf}
lw  <- spdep::nb2listw(spdep::cell2nb(7, 7))
x   <- as.numeric(seq_len(49))
y   <- rev(x)
res <- local_moran_bv(x, y, lw, nsim = 99L, n.cores = 1L, moments = TRUE)
head(res)
attr(res, "quadr")
\dontshow{\}) # examplesIf}
}

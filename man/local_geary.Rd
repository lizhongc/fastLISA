\name{local_geary}
\alias{local_geary}
\title{Univariate Local Geary's C}
\usage{
local_geary(
  x,
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
\item{x}{Numeric vector of length \code{n}.}

\item{listw}{A \code{listw} object from \pkg{spdep}.}

\item{nsim}{Integer; number of permutations. Default \code{999L}. Must be at least 1.}

\item{scale}{Logical; if \code{TRUE} (default), standardise data in R.}

\item{iseed}{Integer seed for RNG, or \code{NULL}.}

\item{p.value}{Numeric significance cutoff. Default \code{0.05}.}

\item{n.cores}{Integer; number of OpenMP threads. Default \code{1L};
set higher to use multiple cores. Ignored on platforms without OpenMP.}

\item{moments}{Logical; if \code{TRUE}, append permutation-distribution
moments \code{E.Ci}, \code{Var.Ci}, \code{Skew.Ci}, and \code{Kurt.Ci}.
Default \code{FALSE}.}

\item{p.method}{Character; permutation pseudo p-value method, \code{"count"}
(default) or \code{"rank"}. \code{"count"} counts permutations at least as
extreme as the observed value; \code{"rank"} uses \pkg{spdep}'s ties-averaged
rank of the observed value. Both return a folded (smaller-tail) value.}
}
\value{
A numeric matrix of class \code{c("localC", "matrix", "array")}
  with columns \code{Ci}, \code{Z.Ci}, and \code{Pr Sim}. When
  \code{moments = TRUE}, the permutation-moment columns are appended.
  It has the following attributes:
  \describe{
    \item{cluster}{A significance-filtered factor with levels Not significant,
      High-High, Low-Low, Other Positive, Negative, Undefined, and Isolated.}
    \item{call}{The matched call.}
  }
}
\description{
\code{local_geary} computes the univariate local Geary's \eqn{C_i}, a
squared-difference Local Indicator of Spatial Association that measures how much a
unit differs from its neighbours. On the sample (\eqn{n-1}) standardised variable
\eqn{z} (when \code{scale = TRUE}) with row-standardised weights \eqn{w^*_{ij}},
\deqn{C_i = \sum_j w^*_{ij} (z_i - z_j)^2 = z_i^2 - 2 z_i\, \mathrm{lag}(z)_i + \mathrm{lag}(z^2)_i.}
A small \eqn{C_i} means \eqn{i} resembles its neighbours (positive spatial
association); a large \eqn{C_i} means it differs from them (negative association,
a spatial outlier).
}
\details{
Inference uses a one-tailed conditional permutation test (\code{nsim} reps): the
observed \eqn{C_i} is compared with the permutation mean to select the tail, and
\deqn{p_i = \frac{t + 1}{\mathrm{nsim} + 1},}
where \eqn{t} counts permuted statistics in that tail. The standardised score is
\eqn{Z.C_i = (C_i - E_{\mathrm{perm}}) / \sqrt{\mathrm{Var}_{\mathrm{perm}}}};
\code{Skew.Ci}/\code{Kurt.Ci} (when \code{moments = TRUE}) follow the \pkg{e1071}
type-3 convention. The \code{cluster} factor splits significant positive
association into \code{High-High}/\code{Low-Low}/\code{Other Positive} and labels
significant dissimilarity \code{Negative}.

Observations with a missing \code{x} value are labelled \code{Undefined} and
observations with no neighbours are labelled \code{Isolated}; both receive
\code{NA} for the p-value, Z-score and moments. The C backend re-seeds its
random number generator per observation, so results are identical for any
\code{n.cores}; \code{n.cores} is ignored when the package is built without
OpenMP.
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
res <- local_geary(x, lw, nsim = 99L, n.cores = 1L)
head(res)
\dontshow{\}) # examplesIf}
}

\name{lmer2lm}
\alias{lmer2lm}

\title{
Coerce objects of class lmerMod or lme to lm
}
\description{Coerce fitted mixed effect models of class lmerMod or lme to an lm-like object using a pseudo data approach (Ritz et al., 2017)
}
\usage{
lmer2lm(modelObject)
}

\arguments{
  \item{modelObject}{a fitted object of class "lmerMod" or "lme", as can be obtained from packages "lme4" and "nlme"}
  }

\details{
  The returned object resembles a fitted model of class "lm" such that its asymptotic representation corresponds to that of the originally fitted mixed model, see Ritz et al. (2017) for details.
}

\references{
    Ritz, C., Pilmann Laursen, R., & Trab Damsgaard, C. (2017). 
  Simultaneous inference for multilevel linear mixed models-with an application to a large-scale school meal study.    \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} \bold{66(2)}, 295--311.
}

\author{Christian Ritz}

%\note{}

%\seealso{}
\examples{
# Variance-covariance matrix using lmer2lm
mSys <- lme(SystolicBP ~ Visit + Gender, random = ~1|ID, data = HealthPromotion, na.action=na.omit)
mDia <- lme(DiastolicBP ~ Visit + Gender, random = ~1|ID, data = HealthPromotion, na.action=na.omit)

mSyslm<-lmer2lm(mSys)
mDialm<-lmer2lm(mDia)

# Variance-covariance matrix for the fixed effects for Visit from both models
m <- mjust(list(mSyslm,mDialm),list("Visit","Visit"))
m$covar
}


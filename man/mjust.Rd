\name{mjust}
\alias{mjust}

\title{
Variance-covariance matrix from multiple marginal models
}
\description{Estimation of the joint variance-covariance matrix from multiple marginal models 
}
\usage{
mjust(modelList, expressions, dataused, level=0.95, seType="san")
}

\arguments{
  \item{modelList}{List of objects of class \code{lm, coxph, nls, drc, lmerMod}, or \code{lme} from which parameters are to be taken}
  \item{expressions}{List of parameters or derived parameters from each model in modelList}
  \item{dataused}{Data set used for all models. The argument is only necessary when modelList includes a model object of type \code{lmerMod} }
  \item{level}{Level of confidence intervals used. Defaults to 0.95}
  \item{seType}{Type of standard errors to report. Default is sandwich variance estimates. Model based estimates may also be chosen using "mod"}
}

\details{
  The function estimates a joint sandwich variance-covariance matrix for all parameters/ derived parameters provided in the expression argument. for details see Pipper et al. 2012 and Jensen et al. 2015.
  Notice: Model objects of class \code{lmerMod} and \code{lme} can have random intercepts and slopes but no nested random effects. For models with several layers of random effects see Ritz et al. 2017. 
  Model objects of class \code{lme} needs to run on complete data sets. Alternatively, use \code{lmer2lm()}.
}

\value{
  A matrix with two columns, one containing BMD and the other containing BMDL.
}

\references{
  
  Jensen, S. M., & Ritz, C. (2015). 
  Simultaneous inference for model averaging of derived parameters. 
  \emph{Risk Analysis} \bold{35(1)}, 68--76.
  
  Jensen, S. M., Pipper, C. B., & Ritz, C. (2015). 
  Evaluation of multi-outcome longitudinal studies. 
  \emph{Statistics in medicine} \bold{34(12)}, 1993--2003.
  
  Pipper, C. B., Ritz, C., & Bisgaard, H. (2012). 
  A versatile method for confirmatory evaluation of the effects of a covariate in multiple models. 
  \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} \bold{61(2)}, 315--326.
  
  Ritz, C., Pilmann Laursen, R., & Trab Damsgaard, C. (2017). 
  Simultaneous inference for multilevel linear mixed models-with an application to a large-scale school meal study.    \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} \bold{66(2)}, 295--311.
}

\author{Signe M. Jensen}

%\note{}

%\seealso{}

\examples{
# Example linear mixed models
require(lme4, quietly=TRUE)

mSys <- lmer(SystolicBP ~ Visit + Gender+(1|ID), data = HealthPromotion)
mDia <- lmer(DiastolicBP ~ Visit + Gender+ (1|ID), data = HealthPromotion)

# Variance-covariance matrix for the fixed effects for Visit from both models
m <- mjust(list(mSys,mDia),list("Visit","Visit"), dataused = HealthPromotion)
m$covar

# Simultaneous confidence intervals
require(multcomp, quietly=TRUE)

beta <- m$coef[,"Estimate"]
Sigma <- m$covar
confint(glht(model = parm(beta, Sigma),alternative="two.sided"))$confint



}
\keyword{variance-covariance matrix}
\keyword{linear mixed model}
\keyword{derived paramters}
\keyword{sandwich variance estimates}

\concept{multiple marginal models}
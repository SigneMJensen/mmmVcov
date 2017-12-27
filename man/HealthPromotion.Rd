\name{HealthPromotion}
\alias{HealthPromotion}
\docType{data}
\title{Health related measurements for university employees enrolled in a workplace health promotion program}
\description{Data contains information on health related parameters for 1365 employees from a public university enrolled in a health promotion program at the workplace. Employees were followed for up till 6 years. The main focus was change in blood pressure over time. For more information see Eng et al. (2016).}
\usage{data(HealthPromotion)}
\format{
  A data frame with 4605 observations on the following 17 variables.
  \describe{
    \item{\code{ID}}{Subject identifier for the 1365 subjects}
    \item{\code{Visit}}{Visit number: 1,...,6}
    \item{\code{Gender}}{Reported as "MALE" or "FEMALE"}
    \item{\code{Ethnicity}}{Reported as "CHINESE", "INDIAN", "MALEY" or "OTHER"}
    \item{\code{Age}}{Age at baseline}
    \item{\code{Occupation}}{Occupation at baseline. Reported as "ACADEMIC" (lecturers, professors and post-doctoral staff), "SUPPORT I" (executives and officers) or "SUPPORT II" (administrative staff, technicians and general workers)}
    \item{\code{Smoking}}{Smoking category at baseline. Reported as "NON-SMOKER","EX-SMOKER" or "SMOKER"}
    \item{\code{BMIGroup}}{BMI group at baseline. Reported as "NORMAL", "OVERWEIGHT" or "OBESE"}
    \item{\code{HypertensionHist}}{History of hypertension at baseline. "YES" or "NO"}
    \item{\code{HypertensionMed0}}{Hypertension medicine at baseline. "YES" or "NO"}
    \item{\code{HypertensionGroup}}{Hypertension group at baseline. "HEALTHY", "AT RISK", "HYPERTENSIVE"}
    \item{\code{Height}}{Height (m)}
    \item{\code{Weight}}{Weight (kg)}
    \item{\code{BMI}}{BMI (kg/m^2)}
    \item{\code{HypertensionMed}}{Hypertension medicine at visit}
    \item{\code{SystolicBP}}{Systolic blod pressure (mmHg)}
    \item{\code{DiastolicBP}}{Diastolic blod pressure (mmHg)}
  }
}
\references{
Eng, J. Y., Moy, F. M., & Bulgiba, A. (2016). Impact of a workplace health promotion program on employees' blood pressure in a public university. PloS one, 11(2), e0148307.
}



\keyword{datasets}

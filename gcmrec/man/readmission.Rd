\name{readmission}
\alias{readmission}
\non_function{}
\title{Rehospitalization colorectal cancer}
\usage{data(readmission)}
\description{
 This contains rehospitalization times after surgery in patients diagnosed with colorectal cancer
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{identificator of each subject. Repeated for each recurrence}
    \item{enum}{which readmission}
    \item{t.start}{start of interval (0 or previous recurrence time)}
    \item{t.stop}{recurrence or censoring time}
    \item{time}{interocurrence or censoring time}
    \item{event}{censoring status. All event are 1 for each subject 
                 excepting last one that it is 0}
    \item{chemo}{Did patient receive chemotherapy? 1: No; 2:Yes}
    \item{sex}{gender: 1:Males 2:Females}
    \item{dukes}{Dukes' tumoral stage: 1:A-B; 2:C 3:D}
    \item{charlson}{Comorbidity Charlson's index. Time-dependent covariate. 0: Index 0; 1: Index 1-2; 3: Index >=3 }
  }
 
}
\source{
 Gonz\'alez, JR., Fernandez, E., Moreno, V. et al. Gender differences in hospital 
 readmission among colorectal cancer patients. Journal of Epidemiology and Community
 Health. In press, 2005.
}

\keyword{datasets}
\keyword{survival}
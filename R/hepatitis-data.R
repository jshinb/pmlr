#' Post-transfusion hepatitis: impact of non-A, non-B hepatitis surrogate tests
#'
#' Data from a randomized double-blind trial to assess whether withholding donor
#' blood positive for the non-A, non-B (NANB) surrogate markers would reduce the
#' frequency of post-transfusion hepatitis.  The dataset contains 4,588 subjects
#' enrolled from 1988 to 1992 into two study groups that received allogenic blood
#' from which units positive for NANB surrogate markers were withheld (n=2,311) or
#' not withheld (n=2,277).  Subjects were followed up for 6 months and assessed
#' for the presence of post-transfusion hepatitis.
#'
#' @usage data(hepatitis)
#'
#' @docType data
#'
#' @format A data frame with 218 rows and the following 6 columns:
#'
#' \tabular{lcl}{
#' \code{city} \tab\tab Subjects were recruited from 3 Canadian Red Cross Society Blood Centres
#' and 13 university-affiliated hopsitals in 3 cities: Toronto, Hamilton and Winnipeg.\cr
#' \code{group} \tab\tab Eligible subjects were assigned to one of two allogenic blood recipient groups.
#' One group received products that had only routine Canadian transfusion-transmissible
#' disease marker screening (no-withhold).
#' The other group received only products that were not positive for NANB surrogate markers (withhold). \cr
#' \code{time} \tab\tab Hepatitis C (HCV) screening was introduced in Canada in May, 1990.
#' Subjects were recruited into the study before (pre) and after (post) the introduction of
#' anti-HCV testing. \cr
#' \code{HCV} \tab\tab Post-transfusion HCV hepatitis present (1) or absent (0). \cr
#' \code{nonABC} \tab\tab Post-transfusion non-A, non-B, non-C hepatitis present (1) or absent (0). \cr
#' \code{counts} \tab\tab Number of subjects.
#'     }
#' @source Blajchman M. A., Bull, S. B. and Feinman S. V. for the Canadian
#' Post-Transfusion Hepatitis Prevention Study Group (1995) Post-transfusion hepatitis:
#' impact of non-A, non-B hepatitis surrogate tests. \emph{The Lancet}, \bold{345}, 21--25.
#' @keywords datasets
"hepatitis"

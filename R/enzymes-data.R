#' Liver Enzyme Data
#'
#' Liver enzyme data collected from 218 patients with liver disease (Plomteux, 1980).
#' The laboratory profile consists of enzymatic activity measured for four liver enzymes:
#' aspartate aminotransferase (AST), alanine aminotransferase (ALT), glutamate dehydrogenase (GLDH)
#' and ornithine carbamyltransferase (OCT).
#'
#' @usage data(enzymes)
#'
#' @docType data
#'
#' @format A data frame with 218 rows and the following 6 columns:
#' \tabular{lcl}{
#' Patient \tab\tab Patient ID\cr
#' Group \tab\tab Four diagnostic groups were considered: acute viral hepatitis (1), persistent chronic hepatitis (2), aggressive chronic hepatitis (3) and post-necrotic cirrhosis (4). \cr
#' AST   \tab\tab Aspartate aminotransferase (U/L) \cr
#' ALT   \tab\tab Alanine aminotransferase (U/L) \cr
#' GLDH  \tab\tab Glutamate dehydrogenase (U/L)\cr
#' OCT   \tab\tab Ornithine carbamyltransferase (U/L)
#'     }
#' @source Albert, A. and Harris, E. K. (1984) \emph{Multivariate Interpretation of Clinical Laboratory Data}. Dekker: New York, Chapter 5, Appendix I.
#' @references Plomteux, G. (1980) Multivariate analysis of an enzyme profile for the differential diagnosis of viral hepatitis. \emph{Clinical Chemistry}, \bold{26}, 1897--1899.
#' @keywords datasets
"enzymes"

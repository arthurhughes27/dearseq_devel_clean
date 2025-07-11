#'
#'@details
#'Analysis of RNA-seq data with variance component score
#'test accounting for data heteroscedasticity through precision weights.
#'Performs gene-wise analysis as well as gene set analysis, including for
#'complex experimental designs such as longitudinal data.
#'
#'\tabular{ll}{
#'Package: \tab dearseq\cr
#'Type: \tab Package\cr
#'Version: \tab 1.13.6\cr
#'Date: \tab 2024-03-06\cr
#'License:\tab \href{http://www.gnu.org/licenses/gpl-2.0.txt}{GPL-2}\cr
#'}
#'The two main functions of the \code{dearseq} package are
#'\code{\link{dear_seq}} and \code{\link{dgsa_seq}}.
#'
#'
#'@references Agniel D & Hejblum BP (2017). Variance component score test for
#'time-course gene set analysis of longitudinal RNA-seq data,
#'\emph{Biostatistics}, 18(4):589-604.
#'\href{https://doi.org/10.1093/biostatistics/kxx005}{DOI: 10.1093/biostatistics/kxx005}.
#'\href{https://arxiv.org/abs/1605.02351}{arXiv:1605.02351}.
#'
#'Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2020). dearseq: a variance
#'component score test for RNA-Seq differential analysis that effectively
#'controls the false discovery rate, \emph{NAR Genomics and Bioinformatics},
#'2(4):lqaa093.
#'\href{https://doi.org/10.1093/nargab/lqaa093}{DOI: 10.1093/nargab/lqaa093}.
#'\href{https://www.biorxiv.org/content/10.1101/635714}{DOI: 10.1101/635714}
#'
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
#' @useDynLib dearseq, .registration = TRUE
## usethis namespace: end
"_PACKAGE"

#' @section Base methods:
#'
#' \describe{
#'   \item{\code{new()}}{Create new instance of this class. Can be called implicitly as \code{ThisSubclassName()}.}
#'
#'   \item{\code{read_eth_files(data_path, labels_path)}}{Read the ETH format files (data and labels).}
#'   \item{\code{read_plink_files(base_path)}}{Read the PLINK format files (with the same \code{base_path} base name).}
#'
#'   \item{\code{write_eth_files(data_path, labels_path)}}{Write the ETH format files (data and labels).}
#'
#'   \item{\code{set_alpha(alpha)}}{Set the significance threshold (FWER).}
#'   \item{\code{get_alpha(alpha)}}{Get the significance threshold (FWER; \code{NULL} if not set).}
#'   \item{\code{set_lmax(lmax)}}{Set the maximum number of SNPs to consider (\code{0} means consider all).}
#'   \item{\code{get_lmax(lmax)}}{Get the maximum number of SNPs to consider (\code{0} means consider all; ; \code{NULL} if not set).}
#'
#'   \item{\code{execute()}}{Execute the search.}
#'
#'   \item{\code{get_result()}}{Return all results of the search.}
#'
#'   \item{\code{write_summary(file_path)}}{Write search summary to a file.}
#'   \item{\code{write_profile(file_path)}}{Write profiling summary to a file.}
#' }

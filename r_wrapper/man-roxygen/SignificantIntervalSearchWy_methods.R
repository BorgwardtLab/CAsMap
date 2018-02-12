#' @section Power-full (WY) methods:
#'
#' \describe{
#'   \item{\code{set_seed(seed)}}{Set random seed (for reproducibility).}
#'
#'   \item{\code{set_n_perm(n_perm)}}{Set number of permutations of the observations used for estimation of intermediate FWERs.}
#'   \item{\code{get_n_perm()}}{Get current number of permutations (\code{NULL} if not set).}
#'   \item{\code{set_perm_file_name()}}{Set file name to write the permutations to (during \code{execute} call).}
#' }

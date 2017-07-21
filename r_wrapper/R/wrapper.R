library(methods)

#new wrapper to go at end

.SignificantFeaturesSearch <- setRefClass(".SignificantFeaturesSearch",
    fields = c('.inst', '.alpha', '.lmax', '.result_available', '.file_loaded'),
    methods = list(

        initialize = function(set_defaults=TRUE)
        {
            .self$.inst <-.self$.create_instance()
            .self$.alpha <- NULL
            .self$.lmax <- NULL
            .self$.result_available <- FALSE
            .self$.file_loaded <- FALSE
            if (set_defaults) {
                .self$set_alpha(0.05)
                .self$set_lmax(0)
            }
        },

        .create_instance = function()
        {
            NULL
        },

        .delete_instance = function(inst)
        {
            NULL
        },

        finalize = function()
        {
            .self$.delete_instance(.self$.inst)
        },

        .check_if_read_is_allowed = function()
        {
        },

        .mark_read_done = function() {
            .self$.result_available <- FALSE
            .self$.file_loaded <- TRUE
        },

        read_eth_files = function(data_path, labels_path, ...)
        {
            .self$.check_if_read_is_allowed()
            .self$.do_read_eth_files(data_path, labels_path, ...)
            .self$.mark_read_done()
        },

        # TODO: split for _int and _iset due to multi-inheritance
        .do_read_eth_files = function(data_path, labels_path, ...) {
            lib_read_eth_files(.self$.inst, data_path, labels_path)
        },

        read_plink_files = function(base_path, ...)
        {
            .self$.check_if_read_is_allowed()
            .self$.do_read_plink_files(base_path, ...)
            .self$.mark_read_done()
        },

        # TODO: split for _int and _iset due to multi-inheritance
        .do_read_plink_files = function(base_path, ...) {
            lib_read_plink_files(.self$.inst, base_path)
        },

        .check_if_write_is_allowed = function() {
            .self$.check_if_files_are_loaded()
        },

        write_eth_files = function(x_file, y_file, ...)
        {
            .self$.check_if_write_is_allowed()
            .self$.write_eth_files(x_file, y_file, ...)
        },

        .check_if_alpha_value_is_allowed = function(x)
        {
            if (!(x>=0 && x<=1)) {
                stop("you need to set alpha to a value between 0 and 1")
            }
        },

        set_alpha = function(alpha)
        {
            .self$.check_if_alpha_value_is_allowed(alpha)
            .self$.alpha <- alpha
            .self$.result_available <- FALSE
        },

        get_alpha = function() {
            return(.self$.alpha)
        },

        .check_if_lmax_value_is_allowed = function(x)
        {
            if (!(x%%1==0 && x>=0)) {
                stop("you need to set lmax to a non-negative integer value")
            }
        },

        set_lmax = function(lmax)
        {
            .self$.check_if_lmax_value_is_allowed(lmax)
            .self$.lmax <- lmax
            .self$.result_available <- FALSE
        },

        get_lmax = function() {
            return(.self$.lmax)
        },

        execute = function()
        {
            .self$.check_if_execute_is_allowed()
            .self$.result_available <- FALSE
            .self$.execute()
            .self$.result_available <- TRUE
        },

        .check_if_result_available = function()
        {
            if (!.self$.result_available) {
                stop("you need to call the execute method first")
            }
        },

        get_result = function()
        {
            .self$.check_if_result_available()
            result <- .self$.get_result()
            result["significance.level"] <-.self$.alpha
            return(result)
        },

        write_summary = function(file_path)
        {
            .self$.check_if_result_available()
            .self$.write_summary(file_path)
        },

        .check_if_files_are_loaded = function()
        {
            if (!.self$.file_loaded) {
                stop("you need to call the read_eth_files / read_plink_files method first")
            }
        },

        .check_if_alpha_is_set = function()
        {
            if (is.null(.self$.alpha)) {
                stop("you need to call the set_alpha method first")
            }
        },

        .check_if_lmax_is_set = function()
        {
            if (is.null(.self$.lmax)) {
                stop("you need to call the set_lmax method first")
            }
        },

        .check_if_execute_is_allowed = function()
        {
            .self$.check_if_alpha_is_set()
            .self$.check_if_lmax_is_set()
            .check_if_files_are_loaded()
        },

        write_profile = function(file_path)
        {
            .self$.check_if_result_available()
            lib_profiler_write_to_file(.self$.inst, file_path)
        }

    )
)

.SignificantIntervalSearch <- setRefClass(".SignificantIntervalSearch",
    contains = c(".SignificantFeaturesSearch"),
    methods = list(

        # .do_read_eth_files = function(data_path, labels_path, ...)
        # {
        #     lib_read_eth_files_int(.self$.inst, data_path, labels_path)
        # },

        # .do_read_plink_files = function(base_path, ...)
        # {
        #     lib_read_plink_files_int(.self$.inst, base_path)
        # },

        .execute = function() {
            lib_execute_int(.self$.inst, .self$.alpha, .self$.lmax)
        },

        .write_eth_files = function(x_file, y_file, ...)
        {
            lib_write_eth_files_int(.self$.inst, x_file, y_file)
        },

        write_filtered_intervals = function(file_path)
        {
            .self$.check_if_result_available()
            lib_filter_intervals_write_to_file(.self$.inst, file_path)
        },

        write_pvals_testable_intervals = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_testable_ints_write_to_file(.self$.inst, file_path)
        },

        write_pvals_significant_intervals = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_significant_ints_write_to_file(.self$.inst, file_path)
        },

        get_significant_intervals = function()
        {
            .self$.check_if_result_available()
            return(lib_get_significant_intervals(.self$.inst))
        },

        get_filtered_intervals = function()
        {
            .self$.check_if_result_available()
            return(lib_get_filtered_intervals(.self$.inst))
        },

        .get_result = function()
        {
            return(lib_get_result_fais(.self$.inst))
        }
    )

)

.SignificantIntervalSearchFais <- setRefClass(".SignificantIntervalSearchFais",
    contains = c(".SignificantIntervalSearch"),
    methods = list(

        .write_summary = function(file_path)
        {
            lib_summary_write_to_file_fais(.self$.inst, file_path)
        }

    )

)

#' Exact fast significant interval search
#'
#' Class for exact significant intervals search with Taron correction for
#' bounding intermediate FWERs.
#'
#' @format \code{\link{RefClass}} object.
#'
#' @template SignificantFeaturesSearch_methods
#' @template SignificantIntervalSearch_methods
SignificantIntervalSearchExact <- setRefClass("SignificantIntervalSearchExact",
    contains = ".SignificantIntervalSearchFais",
    methods = list(
        .create_instance = lib_new_search_e,
        .delete_instance = lib_delete_search_e
    )
)

#' Approximate fast significant interval search
#'
#' Class for approximate significant intervals search with Taron correction for
#' bounding intermediate FWERs.
#'
#' @format \code{\link{RefClass}} object.
#'
#' @template SignificantFeaturesSearch_methods
#' @template SignificantIntervalSearch_methods
SignificantIntervalSearchChi <- setRefClass("SignificantIntervalSearchChi",
    contains = c(".SignificantIntervalSearchFais"),
    methods = list(
        .create_instance = lib_new_search_chi,
        .delete_instance = lib_delete_search_chi
    )
)


.SignificantIntervalSearchWy <- setRefClass(".SignificantIntervalSearchWy",

    contains = c(".SignificantIntervalSearchFais"),
    fields = c(".n_perm"),
    methods = list(

        initialize = function(set_defaults=TRUE, ...) {
            callSuper(set_defaults=set_defaults, ...)
            .self$.n_perm <- NULL
            if (set_defaults) {
                .self$set_n_perm(50)
            }
        },

        .write_summary = function(file_path)
        {
            lib_summary_write_to_file_wy(.self$.inst, file_path)
        },

        .check_if_n_perm_value_is_allowed = function(x)
        {
            if (!(x%%1==0 && x>0)) {
                stop("you need to set n_perm to a positive integer value")
            }
        },

        set_n_perm = function(n_perm) {
            .self$.check_if_n_perm_value_is_allowed(n_perm)
            lib_set_n_perm_wy(.self$.inst, n_perm)
            .self$.n_perm <- n_perm
            .self$.result_available <- FALSE
        },

        get_n_perm = function() {
            return(.self$.n_perm)
        },

        set_seed = function(seed) {
            lib_set_seed_wy(.self$.inst, seed)
            .self$.result_available <- FALSE
        },

        set_perm_file_name = function(file_name) {
            lib_set_perm_file_name_wy(.self$.inst, file_name)
        }
    )
)

#' Exact power-full significant interval search
#'
#' Class for exact significant intervals search with Westfall-Young
#' permutation-based estimation of intermediate FWERs.
#'
#' @format \code{\link{RefClass}} object.
#'
#' @template SignificantFeaturesSearch_methods
#' @template SignificantIntervalSearch_methods
#' @template SignificantIntervalSearchWy_methods
SignificantIntervalSearchWyExact <- setRefClass("SignificantIntervalSearchWyExact",
    contains = c(".SignificantIntervalSearchWy"),
    methods = list(
        .create_instance = lib_new_search_wy_e,
        .delete_instance = lib_delete_search_wy_e
    )
)

#' Approximate power-full significant interval search
#'
#' Class for approximate significant intervals search with Westfall-Young
#' permutation-based estimation of intermediate FWERs.
#'
#' @format \code{\link{RefClass}} object.
#'
#' @template SignificantFeaturesSearch_methods
#' @template SignificantIntervalSearch_methods
#' @template SignificantIntervalSearchWy_methods
SignificantIntervalSearchWyChi <- setRefClass("SignificantIntervalSearchWyChi",
    contains = c(".SignificantIntervalSearchWy"),
    methods = list(
        .create_instance = lib_new_search_wy_chi,
        .delete_instance = lib_delete_search_wy_chi
    )
)



.SignificantFeaturesSearchWithCovariates <- setRefClass(".SignificantFeaturesSearchWithCovariates",
    contains = c(".SignificantFeaturesSearch"),
    fields = c(".cov_loaded"),
    methods = list(

        initialize = function(...) {
            callSuper(...)
            .self$.cov_loaded <- FALSE
        },

        .do_read_eth_files = function(data_path, labels_path, cov_path=NULL, ...) {
            if (is.null(cov_path)) {
                callSuper(data_path, labels_path, ...)
                .self$.cov_loaded <- FALSE
            } else {
                .self$.read_eth_files_with_cov(data_path, labels_path, cov_path)
                .self$.cov_loaded <- TRUE
            }
        },

        .do_read_plink_files = function(base_path, cov_path=NULL, ...) {
            if (is.null(cov_path)) {
                callSuper(base_path, ...)
                .self$.cov_loaded <- FALSE
            } else {
                .self$.read_plink_files_with_cov(base_path, cov_path)
                .self$.cov_loaded <- TRUE
            }
        },

        update_covariates_file = function(cov_path)
        {
            .self$.check_if_files_are_loaded()
            .self$.update_covariates_file(cov_path)
            .self$.result_available <- FALSE
            .self$.cov_loaded <- TRUE
        },

        .write_eth_files = function(x_file, y_file, cov_path=NULL, ...)
        {
            if (is.null(cov_path)) {
                lib_write_eth_files(.self$.inst, x_file, y_file)
            } else {
                .self$.check_if_covariates_are_loaded()
                .self$.write_eth_files_with_cov(x_file, y_file, cov_path)
            }
        },

        .check_if_covariates_are_loaded = function() {
            if (!.self$.cov_loaded) {
                warning("assuming one covariate for all observations; to change covariates call the update_covariates_file method first")
            }
        },

        .check_if_execute_is_allowed = function()
        {
            callSuper()
            .self$.check_if_covariates_are_loaded()
        },

        .get_result = function()
        {
            return(lib_get_result_int(.self$.inst))
        }

    )
)

#' Fast significant interval search with categorical covariates
#'
#' Class for significant intervals search with Cochran-Mantel-Haenszel test,
#' which accounts for categorical covariates.
#'
#' @format \code{\link{RefClass}} object.
#'
#' @template SignificantFeaturesSearch_methods
#' @template SignificantIntervalSearch_methods
#' @template SignificantFeaturesSearchWithCovariates_methods
SignificantIntervalSearchFastCmh <- setRefClass("SignificantIntervalSearchFastCmh",
    # Beware: order matters for calling overloaded covariates methods
    contains = c(".SignificantFeaturesSearchWithCovariates", ".SignificantIntervalSearch"),
    methods = list(
        .create_instance = lib_new_search_fastcmh,
        .delete_instance = lib_delete_search_fastcmh,

        .read_eth_files_with_cov = function(x_file, y_file, cov_path)
        {
            lib_read_eth_files_with_cov_fastcmh(.self$.inst, x_file, y_file, cov_path)
        },

        .read_plink_files_with_cov = function(base_path, cov_path)
        {
            lib_read_plink_files_with_cov_fastcmh(.self$.inst, base_path, cov_path)
        },

        .write_eth_files_with_cov = function(x_file, y_file, cov_path, ...)
        {
            lib_write_eth_files_with_cov_fastcmh(.self$.inst, x_file, y_file, cov_path)
        },

        .update_covariates_file = function(covariates_path)
        {
            lib_read_covariates_file_fastcmh(.self$.inst, covariates_path)
        },

        .write_summary = function(file_path)
        {
            lib_summary_write_to_file_fastcmh(.self$.inst, file_path)
        }

    )
)



.SignificantItemsetSearch <- setRefClass(".SignificantItemsetSearch",
    contains = c(".SignificantFeaturesSearch"),
    methods = list(

        # .do_read_eth_files = function(data_path, labels_path)
        # {
        #     lib_read_eth_files_iset(.self$.inst, data_path, labels_path)
        # },

        # .do_read_plink_files = function(base_path)
        # {
        #     lib_read_plink_files_iset(.self$.inst, base_path)
        # },

        .execute = function() {
            lib_execute_iset(.self$.inst, .self$.alpha, .self$.lmax)
        },

        .write_eth_files = function(x_file, y_file, ...)
        {
            lib_write_eth_files_iset(.self$.inst, x_file, y_file)
        },

        write_pvals_testable_itemsets = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_testable_isets_write_to_file(.self$.inst, file_path)
        },

        write_pvals_significant_itemsets = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_significant_isets_write_to_file(.self$.inst, file_path)
        },

        get_significant_itemsets = function()
        {
            .self$.check_if_result_available()
            return(lib_get_significant_itemsets(.self$.inst))
        },

        .get_result = function()
        {
            return(lib_get_result_iset(.self$.inst))
        }

    )

)

#' Significant itemsets search with categorical covariates
#'
#' Class for significant itemsets search with Cochran-Mantel-Haenszel test,
#' which accounts for categorical covariates.
#'
#' @format \code{\link{RefClass}} object.
#'
#' @template SignificantFeaturesSearch_methods
#' @template SignificantItemsetSearch_methods
#' @template SignificantFeaturesSearchWithCovariates_methods
SignificantItemsetSearchFacs <- setRefClass("SignificantItemsetSearchFacs",
    # Beware: order matters for calling overloaded covariates methods
    contains = c(".SignificantFeaturesSearchWithCovariates", ".SignificantItemsetSearch"),
    methods = list(
        .create_instance = lib_new_search_facs,
        .delete_instance = lib_delete_search_facs,

        .read_eth_files_with_cov = function(x_file, y_file, cov_path)
        {
            lib_read_eth_files_with_cov_facs(.self$.inst, x_file, y_file, cov_path)
        },

        .read_plink_files_with_cov = function(base_path, cov_path)
        {
            lib_read_plink_files_with_cov_facs(.self$.inst, base_path, cov_path)
        },

        .write_eth_files_with_cov = function(x_file, y_file, cov_path)
        {
            lib_write_eth_files_with_cov_facs(.self$.inst, x_file, y_file, cov_path)
        },

        .update_covariates_file = function(cov_path)
        {
            lib_read_covariates_file_facs(.self$.inst, cov_path)
        },

        .get_result = function()
        {
            return(lib_get_result_facs(.self$.inst))
        },

        .write_summary = function(file_path)
        {
            lib_summary_write_to_file_facs(.self$.inst, file_path)
        }

    )
)


#' A method to check value is numeric and in open interval
#' 
#' Checks if a value is numeric and strictly between two other values.
#' 
#' @param x Value to be checked. Needs to be numeric.
#' 
#' @param lower Lower bound. Default value is \code{0}.
#' 
#' @param upper Upper bound. Default value is \code{1}.
#' 
#' @return If numeric, and  strictly greater than \code{lower} and 
#'         strictly smaller than \code{upper}, then return \code{TRUE}. 
#'         Else return \code{FALSE}.
isInOpenInterval <- function(x, lower=0, upper=1){
    inInterval <- TRUE
    if (is.finite(x)){
        if ( (x <= lower) || (x >= upper) )
            inInterval <- FALSE
    } else {
        inInterval <- FALSE
    }
    return (inInterval)
}


#' Check if a variable is boolean or not
#'
#' Checks if a variable is boolean, if not throws an error, otherwise
#' returns boolean.
#'
#' @param var The variable to be checked (if boolean).
#'
#' @param name The name of the variable to appear in any error message.
#'
#' @return If not boolean (or \code{NA}), throws error.
#'         If \code{NA}, return \code{FALSE}. Otherwise return 
#'         boolean value of \code{var}.
checkIsBoolean <- function(var, name){
    if (is.logical(var)) {
        # if NA, return FALSE
        if (is.na(var))
            return(FALSE)

        # otherwise, just return its value (it must be TRUE/FALSE
        # from above check
        return(var)
    } else {
        message <- paste0("Error: ", name, " is not a boolean.")
        stop(message)
    }
}



#' A constructor for a significant pattern search object
#'
#' Initialises an object for the \code{FAIS}, \code{FastCMH} or \code{FACS}
#' methods. Returns the object.
#'
#' @param method A string specifying the method. Use either \code{'fais'}, 
#'               \code{'fastcmh'} or \code{'facs'}. Default is code{''}.
#'               Another way to set the method is to set the 
#'               \code{use_intervals}, \code{use_combinations} and
#'               \code{use_covariate} boolean flags (see below).
#'
#' @param use_intervals Boolean flag to set whether to use intervals
#'                      or combinations. Default is \code{NA}. 
#'                      At least one of 'use_intervals' or 
#'                      'use_combinations' needs to be set to \code{TRUE}
#'                      or \code{FALSE}.
#'
#' @param use_combinations Boolean flag to set whether to use combinations 
#'                         or intervals Default is \code{NA}.
#'                         At least one of 'use_intervals' or 
#'                         'use_combinations' needs to be set to \code{TRUE}
#'                         or \code{FALSE}.
#'
#' @param use_covariate Boolean flag to set whether or not to use a covariate.
#'                      Default is \code{NA}, which is interpreted as 
#'                      \code{FALSE}.
#'
#' @param alpha Value in range \eqn{(0, 1)} (excluding 0 and 1) which sets
#'              significance threshold. If outside range, error is thrown.
#'              Default value is \code{0.05}.
#'
#' @param max_length A value which, for the interval-based methods \code{FAIS} 
#'                   \code{FastCMH}, if the user so desired, 
#'                   sets the maximum interval length. For 
#'                   example, if \code{max_length=5}, only intervals up to
#'                   length 5 will be considered. Negative values will be 
#'                   interpreted as \code{0},
#'                   which means that there is no restriction on the interval 
#'                   length (all possible intervals will be considered).
#'                   Non-integer values will be rounded down to nearest 
#'                   integer. Default value is \code{0}
#'
#' @export
sigpatsearch <- function(method='',
                         use_intervals=NA, 
                         use_combinations=NA, 
                         use_covariate=FALSE,
                         alpha=0.05,
                         max_length=0
                         ){

    #check max_length
    if (is.finite(max_length)){
        max_length <- floor(max_length)
        if (max_length < 0)
            max_length <- 0
    } else { 
        message <- "'max_length' needs to be either 0 or a positive integer."
        stop(message)
    }

    #check alpha
    if (!isInOpenInterval(alpha)){
        message <- "'alpha' needs to be a value strictly in interval (0, 1)."
        stop(message)
    }


    #TODO: add support for LAMP (use_combinations==TRUE and use_covariate=FALSE)
    #TODO: Add examples
    #TODO: check default case
    if ( (method=='') & (is.na(use_intervals)) & (is.na(use_combinations)) ){
        message <- paste0("Need to set at least one of 'method',",
                          " 'use_intervals' or 'use_combinations'.")
        stop(message)
    }

    if (method != ''){
        #convert to lower case
        method <- tolower(method)

        #check if one of 'fais', 'fastcmh' or 'facs'
        correctName <- FALSE

        #fais
        if (method=='fais'){
            correctName <- TRUE
            use_intervals <- TRUE
            use_combinations <- FALSE
            use_covariate <- FALSE
            sig <- SignificantIntervalSearchExact()
            sig$set_alpha(alpha) 
            sig$set_lmax(max_length)
            return(sig)
        }

        #fastcmh
        if (method=='fastcmh'){
            correctName <- TRUE
            use_intervals <- TRUE
            use_combinations <- FALSE
            use_covariate <- TRUE
            sig <- SignificantIntervalSearchFastCmh()
            sig$set_alpha(alpha) 
            sig$set_lmax(max_length)
            return(sig)
        }
    
        #facs
        if (method=='facs'){
            correctName <- TRUE
            use_intervals <- FALSE
            use_combinations <- TRUE
            use_covariate <- TRUE
            sig <- SignificantItemsetSearchFacs()
            sig$set_alpha(alpha) 
            sig$set_lmax(max_length)
            return(sig)
        }

        if (correctName==FALSE){
            message <- paste0("'method' parameter needs to be one of ", 
                              "'fais', 'fastcmh' or 'facs'. Otherwise, ",
                              "set the 'use_intervals' or 'use_combinations' ",
                              "flags.")
            stop(message)
        }


    }

    #now we look at the case that the name as NOT been set, but the flags
    #use_intervals or use_combinations have been set.

    #FIRST: check how many NAs
    na_intervals <- is.na(use_intervals)
    na_combinations <- is.na(use_combinations)
    trueCount <- na_intervals + na_combinations
    if(trueCount==2){
        message <- paste0("Need to set 'method', or at least one of the ",
                          "boolean flags 'use_intervals' and ",
                          "'use_combinations'.")

        stop(message)
    }

    #NOW check flags are all booleans
    #and convert NAs to FALSE 
    use_intervals <- checkIsBoolean(use_intervals, "use_intervals")
    use_combinations <- checkIsBoolean(use_combinations, "use_combinations")
    use_covariate <- checkIsBoolean(use_covariate, "use_covariate")

    #At this point, either use_intervals or use_combinations must be NOT NA
    #NA values receive lower priority
    if (na_intervals){
        use_intervals <- !(use_combinations)
    }
    
    if (na_combinations){
        use_combinations <- !(use_intervals)
    }

    #FALSE==0
    #TRUE==1
    #If all are FALSE, throw error
    #Don't need use_covariate to be set; NA is interpreted as FALSE

    #do fais
    if ( (use_intervals || !use_combinations) & (!use_covariate)){
        sig <- SignificantIntervalSearchExact()
        sig$set_alpha(alpha) 
        sig$set_lmax(max_length)
        return(sig)
    }

    #do facs/lamp
    if ( (use_intervals || !use_combinations) & (use_covariate)){
        sig <- SignificantIntervalSearchFastCmh()
        sig$set_alpha(alpha) 
        sig$set_lmax(max_length)
        return(sig)
    }

    #do facs/lamp
    if ( (use_combinations || !use_intervals) ){
        sig <- SignificantItemsetSearchFacs()
        sig$set_alpha(alpha) 
        sig$set_lmax(max_length)
        return(sig)
    }

    #should never reach this point
    #if here, return error
    message <- paste0("Error: no correct flags were used.")
    stop(message)
    #alternative: could return FAIS object with warning? 
}



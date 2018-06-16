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

        read_eth_files = function(data_path, labels_path, cov_path=NULL, encoding="dominant")
        {
            .self$.check_if_read_is_allowed()
            .self$.do_read_eth_files(data_path, labels_path, cov_path, encoding)
            .self$.mark_read_done()
        },

        # TODO: split for _int and _iset due to multi-inheritance
        .do_read_eth_files = function(data_path, labels_path, cov_path=NULL, encoding="dominant") {
            lib_read_eth_files(.self$.inst, data_path, labels_path, encoding)
        },

        read_plink_files = function(base_path, cov_path=NULL, encoding="dominant")
        {
            .self$.check_if_read_is_allowed()
            .self$.do_read_plink_files(base_path, cov_path, encoding)
            .self$.mark_read_done()
        },

        # TODO: split for _int and _iset due to multi-inheritance
        .do_read_plink_files = function(base_path, cov_path=NULL, encoding="dominant") {
            lib_read_plink_files(.self$.inst, base_path, encoding)
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
            result["target.fwer"] <-.self$.alpha
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

        #NOTE NEW: print method
        show = function(){
            mytype <- "unknown"
            myclasstype = toString(class(.self)[[1]])
            if (myclasstype=="SignificantIntervalSearchExact")
                mytype = "FAIS"
            if (myclasstype=="SignificantIntervalSearchFastCmh")
                mytype = "FastCMH"
            if (myclasstype=="SignificantItemsetSearchFacs")
                mytype = "FACS"

            #extract alpha and lmxax
            myalpha <- toString(.self$.alpha)
            mylmax <- toString(.self$.lmax)

            #output message to be returned
            message1 <- paste0(mytype, " object with:")
            message2 <- paste0(" * alpha = ", myalpha)
            message3 <- paste0(" * lmax = ", mylmax)
            cat(message1, "\n")
            cat(message2, "\n")
            cat(message3, "\n")
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
#' Class for exact significant intervals search with Tarone correction for
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
#' Class for approximate significant intervals search with Tarone correction for
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


.SignificantFeaturesSearchWithCovariates <- setRefClass(".SignificantFeaturesSearchWithCovariates",
    contains = c(".SignificantFeaturesSearch"),
    fields = c(".cov_loaded"),
    methods = list(

        initialize = function(...) {
            callSuper(...)
            .self$.cov_loaded <- FALSE
        },

        .do_read_eth_files = function(data_path, labels_path, cov_path=NULL, encoding="dominant") {
            if (is.null(cov_path)) {
                callSuper(data_path, labels_path, encoding)
                .self$.cov_loaded <- FALSE
            } else {
                .self$.read_eth_files_with_cov(data_path, labels_path, cov_path, encoding)
                .self$.cov_loaded <- TRUE
            }
        },

        .do_read_plink_files = function(base_path, cov_path=NULL, encoding="dominant") {
            if (is.null(cov_path)) {
                callSuper(base_path, encoding)
                .self$.cov_loaded <- FALSE
            } else {
                .self$.read_plink_files_with_cov(base_path, cov_path, encoding)
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
                #warning("assuming one covariate for all observations; to change covariates call the update_covariates_file method first")
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

        .read_eth_files_with_cov = function(x_file, y_file, cov_path, encoding)
        {
            lib_read_eth_files_with_cov_fastcmh(.self$.inst, x_file, y_file, cov_path, encoding)
        },

        .read_plink_files_with_cov = function(base_path, cov_path, encoding)
        {
            lib_read_plink_files_with_cov_fastcmh(.self$.inst, base_path, cov_path, encoding)
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

        .read_eth_files_with_cov = function(x_file, y_file, cov_path, encoding)
        {
            lib_read_eth_files_with_cov_facs(.self$.inst, x_file, y_file, cov_path, encoding)
        },

        .read_plink_files_with_cov = function(base_path, cov_path, encoding)
        {
            lib_read_plink_files_with_cov_facs(.self$.inst, base_path, cov_path, encoding)
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


CASMAP <- setRefClass("CASMAP",
    fields = c('.mode', '.alpha', '.max_comb_size', '.core', '.use_covariates'),
    methods = list(
        initialize = function(mode, alpha=0.05, max_comb_size=0) {
            .self$.mode <- mode
            .self$.alpha <- alpha
            .self$setMaxCombinationSize(max_comb_size)

            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        getMode = function() {
            return(.self$.mode)
        },

        getTargetFWER = function() {
            return(.self$.alpha)
        },

        getMaxCombinationSize = function() {
            return(.self$.max_comb_size)
        },

        isInitialized = function() {
            return(!is.null(.self$.core))
        },

        .checkInitialized = function() {
            if (is.null(.self$.core)){
                stop("Object not initialized or hyperparameters changed since last execution. Please call method readFiles prior to execute.")
            }
        },

        setMode = function(mode){
            if (!is.element(mode, c('regionGWAS', 'higherOrderEpistasis'))){
                stop("Currently implemented modes: < regionGWAS | higherOrderEpistasis >")
            }
            .self$.mode <- mode
            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        setTargetFWER = function(alpha=0.05) {
            #check alpha
            if (!isInOpenInterval(alpha)){
                stop("Target FWER 'alpha' needs to be a value strictly between 0 and 1.")
            }
            .self$.alpha <- alpha
            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        setMaxCombinationSize = function(max_comb_size=0) {
            if (is.finite(max_comb_size)){
                max_comb_size <- floor(max_comb_size)
                if (.self$.mode == 'higherOrderEpistasis' & max_comb_size > 0){
                    print("The current implementation of higher-order epistasis analyses does not support a limited maximum number of interacting variants. The analysis will be carried out for an unlimited order.")
                    max_comb_size <- 0
                }
                if (max_comb_size < 0){
                    max_comb_size <- 0
                }
            } else {
                stop("Maximum combination size 'max_comb_size' needs to be either 0 (unlimited) or a positive integer.")
            }
            .self$.max_comb_size <- max_comb_size
            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        .createCore = function() {
            if (!is.null(.self$.use_covariates)){
                # Instantiate object of the appropriate depending on options
                if (.self$.mode == 'regionGWAS' & !.self$.use_covariates){
                    .self$.core <- SignificantIntervalSearchChi()
                } else if (.self$.mode == 'regionGWAS' & .self$.use_covariates){
                    .self$.core <- SignificantIntervalSearchFastCmh()
                } else if (.self$.mode == 'higherOrderEpistasis'){
                    .self$.core <- SignificantItemsetSearchFacs()
                }

                # Set parameters of the object
                .self$.core$set_alpha(.self$.alpha)
                .self$.core$set_lmax(.self$.max_comb_size)

            } else {
                .self$.core <- NULL
                .self$.use_covariates <- NULL
            }
        },

        readFiles = function(genotype_file=NULL, phenotype_file=NULL, plink_file_root=NULL, covariate_file=NULL, encoding="dominant") {
            # Check whether user decided to use tab-separated text files (binary_format) or PLINK formatted files (plink_format)
            binary_format <- !is.null(genotype_file) & !is.null(phenotype_file)
            plink_format <- !is.null(plink_file_root)
            # At least one of the two must be two, otherwise raise an error
            if (!(binary_format | plink_format)){
                stop("Either plink_file_root or genotype_file and phenotype_file must be specified as arguments.")
            }
            # Check that encoding type is correct
            if (!is.element(encoding, c('dominant', 'recessive'))){
                stop("Currently implemented encodings: < dominant | recessive >")
            }

            # If an additional covariates file was specified, set the object into "CMH mode"
            .self$.use_covariates <- !is.null(covariate_file)

            # Create appropriate "core" object
            .self$.createCore()

            # Give preference to plink_format over binary_format if, by any reason, a user decides to mess around and
            # specify both
            if (plink_format){
                if(.self$.use_covariates){
                    .self$.core$read_plink_files(plink_file_root, covariate_file, encoding)
                } else {
                    .self$.core$read_plink_files(plink_file_root, encoding)
                }
            } else if(binary_format){
                if(.self$.use_covariates){
                    .self$.core$read_eth_files(data_path=genotype_file, labels_path=phenotype_file, cov_path=covariate_file, encoding=encoding)
                } else {
                    .self$.core$read_eth_files(data_path=genotype_file, labels_path=phenotype_file, encoding=encoding)
                }
            } else{
                stop("Either plink_file_root or genotype_file and phenotype_file must be specified as arguments.")
            }
        },

        execute = function(){
            .self$.checkInitialized()
            .self$.core$execute()
        },

        writeSummary = function(path){
            .self$.checkInitialized()
            .self$.core$write_summary(path)
        },

        writeProfile = function(path){
            .self$.checkInitialized()
            .self$.core$write_profile(path)
        },

        writeSignificantRegions = function(path){
            if (.self$.mode != 'regionGWAS'){
                stop("Method writeSignificantInteractions only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            .self$.core$write_pvals_significant_intervals(path)
        },

        writeSignificantClusterRepresentatives = function(path){
            if (.self$.mode != 'regionGWAS'){
                stop("Method writeSignificantClusterRepresentatives only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            .self$.core$write_filtered_intervals(path)
        },

        writeSignificantInteractions = function(path){
            if (.self$.mode != 'higherOrderEpistasis'){
                stop("Method writeSignificantInteractions only available for higher-order epistasis analyses.")
            }
            .self$.checkInitialized()
            .self$.core$write_pvals_significant_itemsets(path)
        },

        getSummary = function(){
            .self$.checkInitialized()
            return(.self$.core$get_result())
        },

        getSignificantRegions = function(){
            if (.self$.mode != 'regionGWAS'){
                stop("Method getSignificantRegions only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            return(.self$.core$get_significant_intervals())
        },

        getSignificantClusterRepresentatives = function(){
            if (.self$.mode != 'regionGWAS'){
                stop("Method getSignificantClusterRepresentatives only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            return(.self$.core$get_filtered_intervals())
        },

        getSignificantInteractions = function(){
            if (.self$.mode != 'higherOrderEpistasis'){
                stop("Method getSignificantInteractions only available for higher-order epistasis analyses.")
            }
            .self$.checkInitialized()
            return(.self$.core$get_significant_itemsets())
        },

        show = function(){
            cat("CASMAP object with:", "\n")
            cat(paste(" * Mode =", .self$.mode), "\n")
            cat(paste(" * Target FWER =", .self$.alpha), "\n")
            cat(paste(" * Maximum combination size =", .self$.max_comb_size), "\n")
            if (!is.null(.self$.core)){
                cat(" * Input files read", "\n")
                cat(paste(" * Covariate =", .self$.use_covariates), "\n")
            } else{
                cat(" * No input files read", "\n")
            }
        }
    )
)


library(sigpatsearch)
context("sigpatsearch")

remove_output_files <- function(itemsets)
{
        files <- c('output_sigints.csv', 'output_summary.txt',
                   'output_timing.txt');

        if (itemsets) {
            files <- c(files,
                       'p_vals_testable_isets.txt',
                       'p_vals_significant_isets.txt');
        } else {
            files <- c(files, 'output_sigints_filtered.corrected.csv',
                       'p_vals_testable_ints.txt',
                       'p_vals_significant_ints.txt');
        }

        for (file_name in files)
        {
            if (file.exists(file_name))
                file.remove(file_name);
        }

}

check_files <- function(itemsets, rm.files=TRUE)
{
        expect_true(file.exists('output_summary.txt'));
        expect_true(file.exists('output_timing.txt'));
        if (itemsets) {
            expect_true(file.exists('p_vals_testable_isets.txt'));
            expect_true(file.exists('p_vals_significant_isets.txt'));
        } else {
            expect_true(file.exists('output_sigints_filtered.corrected.csv'));
            expect_true(file.exists('p_vals_testable_ints.txt'));
            expect_true(file.exists('p_vals_significant_ints.txt'));
        }
        if (rm.files) {
            remove_output_files(itemsets);
        }
}

test_write_wrapper <- function(inst, summary_write_to_file_fun, itemsets=FALSE) {
    summary_write_to_file_fun(inst, "output_summary.txt");
    lib_profiler_write_to_file(inst, "output_timing.txt");
    if (itemsets) {
        lib_pvals_testable_isets_write_to_file(inst, "p_vals_testable_isets.txt");
        lib_pvals_significant_isets_write_to_file(inst, "p_vals_significant_isets.txt");
    } else {
        lib_filter_intervals_write_to_file(inst, "output_sigints_filtered.corrected.csv");
        lib_pvals_testable_ints_write_to_file(inst, "p_vals_testable_ints.txt");
        lib_pvals_significant_ints_write_to_file(inst, "p_vals_significant_ints.txt");
    }
    check_files(itemsets);
}


.create_check_features_fun <- function(df_colnames) {
    ret <- function(df, n, ...) {
        expect_equal(nrow(df), n, ...);
        expect_true(all(names(df) == df_colnames));
    }
    return(ret)
}
.intervals_colnames = c('start', 'end', 'pvalue')
check_significant_intervals <- .create_check_features_fun(.intervals_colnames)
check_filtered_intervals <- .create_check_features_fun(.intervals_colnames)
.itemsets_colnames = c('itemsets', 'pvalue')
check_significant_itemsets <- .create_check_features_fun(.itemsets_colnames)


test_wrapper <- function(inst, execute_fun, summary_write_to_file_fun, n_expected, itemsets=FALSE, data_file="data.txt", label_file="label.txt")
{
        lib_read_eth_files(inst, data_file, label_file);
        execute_fun(inst, 0.04, 0);
        test_write_wrapper(inst, summary_write_to_file_fun, itemsets=itemsets);

        if (itemsets) {
            df <- lib_get_significant_itemsets(inst);

            check_significant_itemsets(df, n_expected);
        } else {
            df <- lib_get_significant_intervals(inst);
            df2 <- lib_get_filtered_intervals(inst);

            check_significant_intervals(df, n_expected);
            check_filtered_intervals(df2, 2);
        }
}
test_wrapper_fais <- function(inst, ...)
{
    test_wrapper(inst, lib_execute_int, lib_summary_write_to_file_fais, ...)
}
test_that("check wrapper functions exact", {

        inst <- lib_new_search_e();
        test_wrapper_fais(inst, 151)
        lib_delete_search_e(inst);

})
test_that("check wrapper functions chi", {

        inst <- lib_new_search_chi();
        test_wrapper_fais(inst, 158)
        lib_delete_search_chi(inst);

})

test_wrapper_cov <- function(inst, read_covariates_file_fun, execute_fun, summary_write_to_file_fun, n_expected, itemsets=FALSE, data_file="data.txt", label_file="label.txt", cov_file="cov.txt")
{
    lib_read_eth_files(inst, data_file, label_file);
    read_covariates_file_fun(inst, cov_file);
    execute_fun(inst, 0.05, 0);
    test_write_wrapper(inst, summary_write_to_file_fun, itemsets=itemsets);

    if (itemsets) {
        df <- lib_get_significant_itemsets(inst);

        check_significant_itemsets(df, n_expected);
    } else {
        df <- lib_get_significant_intervals(inst);
        df2 <- lib_get_filtered_intervals(inst);

        check_significant_intervals(df, n_expected);
        check_filtered_intervals(df2, 1);
    }
}
test_that("check wrapper functions fastcmh", {

    inst <- lib_new_search_fastcmh();
    data_file="cov_data.txt"; label_file="cov_label.txt";
    test_wrapper(inst, lib_execute_int, lib_summary_write_to_file_fastcmh, 4, data_file=data_file, label_file=label_file);
    test_wrapper_cov(inst, lib_read_covariates_file_fastcmh, lib_execute_int, lib_summary_write_to_file_fastcmh, 5, data_file=data_file, label_file=label_file, cov_file="cov.txt");
    lib_delete_search_fastcmh(inst);
})
test_that("check wrapper functions facs", {

    inst <- lib_new_search_facs();
    data_file="tictactoe_data.txt"; label_file="tictactoe_labels.txt";
    cov_file="tictactoe_covariates.txt";
    test_wrapper(inst, lib_execute_iset, lib_summary_write_to_file_facs, 362, itemsets=TRUE, data_file=data_file, label_file=label_file);
    test_wrapper_cov(inst, lib_read_covariates_file_facs, lib_execute_iset, lib_summary_write_to_file_facs, 365, itemsets=TRUE, data_file=data_file, label_file=label_file, cov_file=cov_file);
    lib_delete_search_facs(inst);
})

test_wrapper_wy <- function(inst, n_expected)
{
        lib_read_eth_files(inst, "data.txt", "label.txt");
        lib_set_n_perm_wy(inst, 4);
        lib_set_seed_wy(inst, 1);
        lib_execute_int(inst, 0.25, 5);
        test_write_wrapper(inst, lib_summary_write_to_file_wy);

        df <- lib_get_significant_intervals(inst);
        df2 <- lib_get_filtered_intervals(inst);

        check_significant_intervals(df, n_expected);
        check_filtered_intervals(df2, 2);
}
test_that("check wrapper functions wy exact", {

        inst <- lib_new_search_wy_e();
        test_wrapper_wy(inst, 60)
        lib_delete_search_wy_e(inst);

})
test_that("check wrapper functions wy chi", {

        inst <- lib_new_search_wy_chi();
        test_wrapper_wy(inst, 60)
        lib_delete_search_wy_chi(inst);

})


test_write_class <- function(search, itemsets=FALSE) {
    search$write_summary("output_summary.txt");
    search$write_profile("output_timing.txt");
    if (itemsets) {
        search$write_pvals_testable_itemsets("p_vals_testable_isets.txt");
        search$write_pvals_significant_itemsets("p_vals_significant_isets.txt");
    } else {
        search$write_filtered_intervals("output_sigints_filtered.corrected.csv");
        search$write_pvals_testable_intervals("p_vals_testable_ints.txt");
        search$write_pvals_significant_intervals("p_vals_significant_ints.txt");
    }
    check_files(itemsets);
}


test_class_impl <- function(search_cls, n_expected, data_file="data.txt", label_file="label.txt") {

    search = search_cls(set_defaults=FALSE)

    # wrong files
    expect_error(search$read_eth_files("aaargh_data.txt", "label.txt"));
    expect_error(search$read_eth_files("data.txt", "aaargh_label.txt"));

    # set_alpha and set_lmax not called yet:
    expect_error(search$execute());

    expect_error(search$set_alpha(-0.01))
    expect_error(search$set_alpha(1.01))
    search$set_alpha(0.04);
    # set_lmax not called yet:
    expect_error(search$execute());
    expect_error(search$write_summary("output_summary.txt"));
    expect_error(search$write_profile("output_timing.txt"));
    expect_error(search$write_filtered_intervals("output_sigints_filtered.corrected.csv"));
    expect_error(search$write_pvals_testable_intervals("p_vals_testable_ints.txt"));
    expect_error(search$write_pvals_significant_intervals("p_vals_significant_ints.txt"));

    expect_error(search$set_lmax(-1))
    expect_error(search$set_lmax(0.1))
    search$set_lmax(0);

    search$read_eth_files(data_file, label_file);
    search$execute();
    test_write_class(search);

    df <- search$get_significant_intervals();
    df2 <- search$get_filtered_intervals();

    check_significant_intervals(df, n_expected)
    check_filtered_intervals(df2, 2)
}



test_that("check class implementation exact", {

    test_class_impl(SignificantIntervalSearchExact, 151);

})

test_that("check class implementation chi", {

    test_class_impl(SignificantIntervalSearchChi, 158);

})


test_class_impl_wy <- function(search_cls, n_expected) {

    search = search_cls(set_defaults=FALSE)

    # without reading files execute does not work:
    expect_error(search$execute());

    search$read_eth_files("data.txt", "label.txt");

    # without setting alpha execute does not work:
    expect_error(search$execute());

    search$set_lmax(5);
    expect_error(search$execute());

    search$set_alpha(0.25);

    search$set_seed(1);

    # without setting n_perm execute does not work:
    expect_error(search$execute());

    expect_error(search$set_n_perm(0));
    expect_error(search$set_n_perm(0.1));
    search$set_n_perm(4);

    # now we have everything initialized so that execute works:
    search$execute();
    test_write_class(search);

    df <- search$get_significant_intervals();
    df2 <- search$get_filtered_intervals();

    check_significant_intervals(df, n_expected)
    check_filtered_intervals(df2, 2)
}


test_that("check class implementation wy exact", {

    test_class_impl_wy(SignificantIntervalSearchWyExact, 60);
})


test_that("check class implementation wy chi", {

    test_class_impl_wy(SignificantIntervalSearchWyChi, 60);
})



test_class_impl_cov <- function(search_cls, n_expected, itemsets=FALSE,
        data_file="cov_data.txt", label_file="cov_label.txt",
        cov_file="cov.txt") {

    search = search_cls(set_defaults=FALSE)

    # without reading files execute does not work:
    expect_error(search$execute());
    search$read_eth_files(data_file, label_file);

    # without setting alpha execute does not work:
    expect_error(search$execute());
    search$set_alpha(0.05);
    expect_error(search$execute());
    search$set_lmax(0);

    # w/o covariates read execute warns about default setting
    expect_warning(search$execute());
    test_write_class(search, itemsets=itemsets);

    if (itemsets) {
        df <- search$get_significant_itemsets();

        check_significant_itemsets(df, n_expected+1);
    } else {
        df <- search$get_significant_intervals();
        df2 <- search$get_filtered_intervals();

        check_significant_intervals(df, n_expected-1)
        check_filtered_intervals(df2, 2)
    }

    # wrong file
    expect_error(search$read_eth_files(data_file, label_file, "aaargh_cov.txt"));
    expect_error(search$update_covariates_file("aaargh_cov.txt"));
    search$update_covariates_file(cov_file);

    search$execute();
    test_write_class(search, itemsets=itemsets);

    if (itemsets) {
        df <- search$get_significant_itemsets();

        check_significant_itemsets(df, n_expected);
    } else {
        df <- search$get_significant_intervals();
        df2 <- search$get_filtered_intervals();

        check_significant_intervals(df, n_expected)
        check_filtered_intervals(df2, 1)
    }
}

test_that("check class implementation fastcmh", {

    data_file="cov_data.txt"; label_file="cov_label.txt"; cov_file="cov.txt"
    test_class_impl_cov(SignificantIntervalSearchFastCmh, 5, data_file=data_file, label_file=label_file,
        cov_file=cov_file)
})

test_that("check class implementation facs", {

    data_file="tictactoe_data.txt"; label_file="tictactoe_labels.txt"
    cov_file="tictactoe_covariates.txt"
    test_class_impl_cov(SignificantItemsetSearchFacs, 365, itemsets=TRUE, data_file=data_file,
        label_file=label_file, cov_file=cov_file)
})

#TODO: Uncomment
#cat("here")

# test_that("check plink reading", {
# 
#     search <- SignificantIntervalSearchWyChi();
# 
#     search$set_n_perm(4);
#     alpha <- 0.25;
#     search$set_alpha(alpha);
#     search$set_lmax(5);
#     search$set_seed(0);
# 
#     search$read_plink_files("sample_data");
#     # now we have everything initialized so that execute works:
#     search$execute();
#     test_write_class(search);
# 
#     cat("there")
# 
#     region.offset = 22;
#     expected <- list(
#         n.sig.int = 293,
#         n.sig.int.clustered = 11,
#         int.processed = 549990,
#         int.testable = 0,
#         testability.threshold = 7.09e-7,
#         significance.level = alpha,
#         corrected.significance.level = 1.98e-6,
#         region = list(
#             start = c(region.offset,500),
#             end = c(500, 1000-region.offset)
#         )
#     );
# 
#     df <- search$get_significant_intervals();
#     df2 <- search$get_filtered_intervals();
# 
#     # OPT add tolerance for randomness
#     check_significant_intervals(df, expected$n.sig.int) #, tolerance=5e-2, scale=expected$n.sig.int)
#     check_filtered_intervals(df2, expected$n.sig.int.clustered) #, tolerance=5e-2, scale=expected$n.sig.int)
# 
#     result <-search$get_result()
# 
#     # here exact check: reported results size vs. search fields size
#     check_significant_intervals(result$sig.int, nrow(df))
#     check_filtered_intervals(result$sig.int.clustered, nrow(df2))
# 
#     expect_equal(result$int.processed, expected$int.processed);
#     expect_equal(result$int.testable, expected$int.testable);
#     expect_equal(result$testability.threshold, expected$testability.threshold, tolerance=1e-3, scale=expected$testability.threshold);
#     expect_equal(result$corrected.significance.level, expected$corrected.significance.level, tolerance=1e-3, scale=expected$corrected.significance.level);
#     expect_equal(result$significance.level, expected$significance.level, tolerance=1e-3, scale=expected$corrected.significance.level);
#     region <- result$testability.region;
#     # OPT add tolerance for randomness: far ends only
#     expect_equal(region$start[1], expected$region$start[1]) #, tolerance=1e-1, scale=region.offset)
#     expect_equal(region$start[2], expected$region$start[2])
#     expect_equal(region$end[1], expected$region$end[1])
#     expect_equal(region$end[2], expected$region$end[2]) #, tolerance=1e-1, scale=region.offset)
# })

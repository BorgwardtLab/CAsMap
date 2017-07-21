library(sigpatsearch)

search <- SignificantIntervalSearchWyChi$new();
search$set_n_perm(4);
search$set_seed(1);
search$set_alpha(.3);
search$set_lmax(5);
search$read_eth_files("data.txt", "label.txt");
# now we have everything initialized so that execute works:
search$execute("output", "");
search$write_summary("output_summary.txt");
search$write_profile("output_timing.txt");
search$write_filtered_intervals("output_sigints_filtered.corrected.csv");
search$write_pvals_testable_intervals("p_vals_testable_ints.txt");
search$write_pvals_significant_intervals("p_vals_significant_ints.txt");
df <- search$get_significant_intervals();
df2 <- search$get_filtered_intervals();

result <- list(iv=df, clustered=df2)

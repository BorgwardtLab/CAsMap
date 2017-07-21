search <- SignificantIntervalSearchWyChi$new();

search$set_n_perm(4);
search$set_alpha(0.25);
search$set_lmax(5);
search$set_seed(5);

# search$read_plink_files("sample_data");
search$read_eth_files("tests/testthat/data.txt", "tests/testthat/label.txt");
# now we have everything initialized so that execute works:
search$execute("output", "");
result <-search$get_result()

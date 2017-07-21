context("Testing new wrapper for sigpatsearch")

library(sigpatsearch)

test_that("Wrapper has 'method=fais'", {
    sig <- sigpatsearch(method='fais')
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchExact"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has 'method=fastcmh'", {
    sig <- sigpatsearch(method='fastcmh')
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchFastCmh"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has 'method=facs'", {
    sig <- sigpatsearch(method='facs')
    sigclass <- class(sig)[[1]]
    solution <- "SignificantItemsetSearchFacs"
    expect_equal(sigclass, solution)
})


test_that("Error when wrong name: 'method=xyz'", {
    expect_error(sig <- sigpatsearch(method='xyz'), "needs to be one")
})


test_that("Error when no parameters set", {
    expect_error(sig <- sigpatsearch(), "Need to set at least one of ")
})


#----------------------------------------------------------------------#
# max_length checks

test_that("Error max_length not finite", {
    expect_error(sig <- sigpatsearch(method='fais', max_length=1/0), 
                 "needs to be either 0 or ")
})


test_that("Error max_length not numeric", {
    expect_error(sig <- sigpatsearch(method='fais', max_length='a'), 
                 "needs to be either 0 or ")
})


test_that("No error when max_length negative", {
    expect_error(sig <- sigpatsearch(method='fais', max_length=-1), NA) 
})


test_that("No error when max_length double", {
    expect_error(sig <- sigpatsearch(method='fais', max_length=2.5), NA) 
})


#----------------------------------------------------------------------#
# alpha checks

test_that("Error when alpha  not numeric", {
    expect_error(sig <- sigpatsearch(method='fais', alpha=NA), 
                 "needs to be a value strictly") 
})


test_that("Error when alpha  not numeric", {
    expect_error(sig <- sigpatsearch(method='fais', alpha='a'), 
                 "needs to be a value strictly") 
})


test_that("Error when alpha outside (0, 1)", {
    expect_error(sig <- sigpatsearch(method='fais', alpha=2), 
                 "needs to be a value strictly") 
})


test_that("Error when alpha outside (0, 1)", {
    expect_error(sig <- sigpatsearch(method='fais', alpha=0), 
                 "needs to be a value strictly") 
})


test_that("Error when alpha outside (0, 1)", {
    expect_error(sig <- sigpatsearch(method='fais', alpha=1), 
                 "needs to be a value strictly") 
})


test_that("No error when alpha outside (0, 1)", {
    expect_error(sig <- sigpatsearch(method='fais', alpha=0.0001), NA) 
})


#----------------------------------------------------------------------#
# checkIsBoolean

test_that("checkIsBoolean works for TRUE", {
    x <- TRUE
    expect_error(b <- checkIsBoolean(x, "x"), NA)
})


test_that("checkIsBoolean works for FALSE", {
    x <- FALSE
    expect_error(b <- checkIsBoolean(x, "x"), NA)
})


test_that("checkIsBoolean throws error for numeric", {
    x <- 1
    expect_error(b <- checkIsBoolean(x, "x"), "not a boolean")
})


test_that("checkIsBoolean throws error for string", {
    x <- "x"
    expect_error(b <- checkIsBoolean(x, "x"), "not a boolean")
})


test_that("checkIsBoolean does NOT throw error NA", {
    x <- NA
    expect_error(b <- checkIsBoolean(x, "x"), NA)
})


#----------------------------------------------------------------------#
#no errors...

test_that("checkIsBoolean returns FALSE for NA", {
    x <- NA
    expect_equal(checkIsBoolean(x, "x"), FALSE)
})


test_that("checkIsBoolean returns FALSE for FALSE", {
    x <- FALSE 
    expect_equal(checkIsBoolean(x, "x"), FALSE)
})


test_that("checkIsBoolean returns TRUE for TRUE", {
    x <- TRUE 
    expect_equal(checkIsBoolean(x, "x"), TRUE)
})

#----------------------------------------------------------------------#
#now using flags

test_that("Wrapper has flags for fais1: 
           use_intervals=TRUE, use_covariate=FALSE", {
    sig <- sigpatsearch(use_intervals=TRUE, use_covariate=FALSE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchExact"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for fais2: 
           use_combinations=FALSE, use_covariate=FALSE", {
    sig <- sigpatsearch(use_combinations=FALSE, use_covariate=FALSE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchExact"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for fais3: 
           use_combinations=FALSE", {
    sig <- sigpatsearch(use_combinations=FALSE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchExact"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for fais4: 
           use_intervals=TRUE", {
    sig <- sigpatsearch(use_intervals=TRUE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchExact"
    expect_equal(sigclass, solution)
})


#----------------------------------------------------------------------#
#let's see what happens if use_intervals and use_combinations are the
#same - it should go to FAIS, the way the code is structured


test_that("Wrapper has flags that should give fais1: 
           use_intervals=TRUE, use_combinations=TRUE", {
    sig <- sigpatsearch(use_intervals=TRUE, use_combinations=TRUE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchExact"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags that should give fais2: 
           use_intervals=FALSE, use_combinations=FALSE", {
    sig <- sigpatsearch(use_intervals=FALSE, use_combinations=FALSE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchExact"
    expect_equal(sigclass, solution)
})


#----------------------------------------------------------------------#
test_that("Wrapper has flags for fastcmh1: 
           use_intervals=TRUE, use_covariate=TRUE", {
    sig <- sigpatsearch(use_intervals=TRUE, use_covariate=TRUE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchFastCmh"
    expect_equal(sigclass, solution)
})



test_that("Wrapper has flags for fastcmh2: 
           use_combinations=FALSE, use_covariate=TRUE", {
    sig <- sigpatsearch(use_combinations=FALSE, use_covariate=TRUE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantIntervalSearchFastCmh"
    expect_equal(sigclass, solution)
})


#----------------------------------------------------------------------#
#FACS

test_that("Wrapper has flags for facs1: 
           use_intervals=FALSE, use_covariate=NA", {
    sig <- sigpatsearch(use_intervals=FALSE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantItemsetSearchFacs"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for facs2: 
           use_intervals=FALSE, use_covariate=FALSE", {
    sig <- sigpatsearch(use_intervals=FALSE, use_covariate=FALSE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantItemsetSearchFacs"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for facs3: 
           use_intervals=FALSE, use_covariate=TRUE", {
    sig <- sigpatsearch(use_intervals=FALSE, use_covariate=TRUE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantItemsetSearchFacs"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for facs4: 
           use_combinations=TRUE, use_covariate=NA", {
    sig <- sigpatsearch(use_combinations=TRUE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantItemsetSearchFacs"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for facs4: 
           use_combinations=TRUE, use_covariate=FALSE", {
    sig <- sigpatsearch(use_combinations=TRUE, use_covariate=FALSE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantItemsetSearchFacs"
    expect_equal(sigclass, solution)
})


test_that("Wrapper has flags for facs5: 
           use_intervals=FALSE, use_covariate=TRUE", {
    sig <- sigpatsearch(use_combinations=TRUE, use_covariate=TRUE)
    sigclass <- class(sig)[[1]]
    solution <- "SignificantItemsetSearchFacs"
    expect_equal(sigclass, solution)
})



#----------------------------------------------------------------------#
#force error message with bad flags


test_that("Wrapper has no correct flags 1:", {
    error_message <- "Need to set"
    expect_error(sig <- sigpatsearch(),
                 error_message)
})


test_that("Wrapper has no correct flags 2:", {
    error_message <- "Need to set"
    expect_error(sig <- sigpatsearch(use_covariate=TRUE),
                 error_message)
})


test_that("Wrapper has no correct flags 3:", {
    error_message <- "is not a boolean"
    expect_error(sig <- sigpatsearch(use_combinations="x", use_covariate=TRUE),
                 error_message)
})


test_that("Wrapper has no correct flags 4:", {
    error_message <- "Need to set"
    expect_error(sig <- sigpatsearch(use_intervals=NA, use_covariate=TRUE),
                 error_message)
})



test_that("Wrapper has no correct flags 5:", {
    error_message <- "is not a boolean"
    expect_error(sig <- sigpatsearch(use_combinations="x", use_covariate=TRUE),
                 error_message)
})


test_that("Wrapper has no correct flags 6:", {
    error_message <- "is not a boolean"
    expect_error(sig <- sigpatsearch(use_combinations=6, use_covariate=TRUE),
                 error_message)
})

import unittest
from sigpatsearch import createSigPatSearch

#checkIsBoolean and isInOpenInterval are private with __
from sigpatsearch import __checkIsBoolean as checkIsBoolean
from sigpatsearch import __isInOpenInterval as isInOpenInterval

TYPE_STRING_FASTCMH = '_SignificantIntervalSearchFastCmh'
TYPE_STRING_FACS = '_SignificantItemsetSearchFacs'

class WrapperTests(unittest.TestCase):


#----------------------------------------------------------------------#
# max_length checks
    def test_maxlength1(self):
        '''testing max_length error - not finite
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, max_length=float('nan') )


    def test_maxlength2(self):
        '''testing max_length error - not finite
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, max_length=float('inf') )


    def test_maxlength3(self):
        '''testing max_length error - not numeric
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, max_length='a')


    def test_maxlength4(self):
        '''testing max_length error - None
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, max_length=None)


    def test_maxlength5(self):
        '''testing max_length - no error when negative
        '''
        try:
            createSigPatSearch(use_intervals=True, max_length=-1)
        except ExceptionType:
            self.fail("createSigPatSearch raised ExceptionType unexpectedly!")


    def test_maxlength6(self):
        '''testing max_length - no error when double
        '''
        try:
            createSigPatSearch(use_intervals=True, max_length=2.5)
        except ExceptionType:
            self.fail("createSigPatSearch raised ExceptionType unexpectedly!")

#----------------------------------------------------------------------#
# alpha checks

    def test_alpha1(self):
        '''testing alpha error - None
        '''
        #print(isInOpenInterval('a'))
        #print(isInOpenInterval(-1))
        #print(isInOpenInterval(float('nan')))

        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, alpha=None)


    def test_alpha2(self):
        '''testing alpha error - not numeric
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, alpha='a')


    def test_alpha3(self):
        '''testing alpha error - outside of interval
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, alpha=2)


    def test_alpha4(self):
        '''testing alpha error - outside of interval (0 is the border)
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, alpha=0)


    def test_alpha5(self):
        '''testing alpha error - outside of interval (1 is the border)
        '''
        with self.assertRaises(ValueError):
            createSigPatSearch(use_intervals=True, alpha=1)


    def test_alpha6(self):
        '''testing alpha no error
        '''
        try:
            createSigPatSearch(use_intervals=True, alpha=0.0001)
        except ExceptionType:
            self.fail("createSigPatSearch raised ExceptionType unexpectedly!")

#----------------------------------------------------------------------#
# checkIsBoolean

    def test_checkIsBoolean1(self):
        '''Check no error for True without error
        '''
        x = True
        try:
            b = checkIsBoolean(x, "x")
        except ExceptionType:
            self.fail("checkIsBoolean raised ExceptionType unexpectedly!")


    def test_checkIsBoolean2(self):
        '''Check no error for False without error
        '''
        x = False
        try:
            b = checkIsBoolean(x, "x")
        except ExceptionType:
            self.fail("checkIsBoolean raised ExceptionType unexpectedly!")


    def test_checkIsBoolean3(self):
        '''Check throws error for numeric which is not boolean
        '''
        x = 1.2
        with self.assertRaises(ValueError):
            b = checkIsBoolean(x, "x")


    def test_checkIsBoolean4(self):
        '''Check no error for numeric value which is boolean 1
        '''
        x = 1
        try:
            b = checkIsBoolean(x, "x")
        except ExceptionType:
            self.fail("checkIsBoolean raised ExceptionType unexpectedly!")


    def test_checkIsBoolean5(self):
        '''Check no error for numeric value which is boolean 0.0
        '''
        x = 0.0
        try:
            b = checkIsBoolean(x, "x")
        except ExceptionType:
            self.fail("checkIsBoolean raised ExceptionType unexpectedly!")


    def test_checkIsBoolean6(self):
        '''Check throws error for string
        '''
        x = "x"
        with self.assertRaises(ValueError):
            b = checkIsBoolean(x, "x")

    def test_checkIsBoolean7(self):
        '''Check no error for None
        '''
        x = None
        try:
            b = checkIsBoolean(x, "x")
        except ExceptionType:
            self.fail("checkIsBoolean raised ExceptionType unexpectedly!")


#----------------------------------------------------------------------#
# checkIsBoolean - no errors, now check return values

    def test_checkIsBoolean8(self):
        '''Check returns False for None
        '''
        x = None
        b = checkIsBoolean(x, "x")
        self.assertEqual(b, False)


    def test_checkIsBoolean9(self):
        '''Check returns False for False
        '''
        x = False
        b = checkIsBoolean(x, "x")
        self.assertEqual(b, False)


    def test_checkIsBoolean9_2(self):
        '''Check returns True for True
        '''
        x = True
        b = checkIsBoolean(x, "x")
        self.assertEqual(b, True)

#----------------------------------------------------------------------#
# Initialising using method

    def test_createSigPatSearch_Method_fais(self):
        '''Testing creation of object where method is specified
           In this case, it is fais
        '''
        sig = createSigPatSearch(method='fAiS')
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchExact"
        self.assertEqual(sigclass, solution)
       

    def test_createSigPatSearch_Method_fastcmh(self):
        '''Testing creation of object where method is specified
           In this case, it is fastcmh
        '''
        sig = createSigPatSearch(method='fastcmh')
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchFastCmh"
        self.assertEqual(sigclass, solution)


    def test_createSigPatSearch_Method_facs(self):
        '''Testing creation of object where method is specified
           In this case, it is fastcmh
        '''
        sig = createSigPatSearch(method='FAcS')
        sigclass = type(sig).__name__
        solution = "_SignificantItemsetSearchFacs"
        self.assertEqual(sigclass, solution)


    def test_createSigPatSearch_Method_error1(self):
        '''Testing creation of object where method is specified
           Throws an error, wrong name
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch(method='xyz')


    def test_createSigPatSearch_Method_error2(self):
        '''Testing creation of object where method is specified
           Throws an error, no parameters set
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch()



#----------------------------------------------------------------------#
# Now using flags

    def test_create_flags_fais1(self):
        '''Wrapper has flags for fais1:
           use_intervals=True, use_covariate=False"
        '''
        sig = createSigPatSearch(use_intervals=True, use_covariate=False)
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchExact"
        self.assertEqual(sigclass, solution)


    def test_create_flags_fais2(self):
        '''Wrapper has flags for fais2:
           use_combinations=False, use_covariate=False"
        '''
        sig = createSigPatSearch(use_combinations=False, use_covariate=False)
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchExact"
        self.assertEqual(sigclass, solution)


    def test_create_flags_fais3(self):
        '''Wrapper has flags for fais3:
           use_intervals=True
        '''
        sig = createSigPatSearch(use_intervals=True)
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchExact"
        self.assertEqual(sigclass, solution)


#----------------------------------------------------------------------#
#let's see what happens if use_intervals and use_combinations are the
#same - it should go to FAIS, the way the code is structured


    def test_create_flags_fais4(self):
        '''Wrapper has flags for fais4:
           use_intervals=True, use_combinations=True
        '''
        sig = createSigPatSearch(use_intervals=True, use_combinations=True)
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchExact"
        self.assertEqual(sigclass, solution)


    def test_create_flags_fais5(self):
        '''Wrapper has flags for fais5:
           use_intervals=False, use_combinations=False
        '''
        sig = createSigPatSearch(use_intervals=False, use_combinations=False)
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchExact"
        self.assertEqual(sigclass, solution)


#----------------------------------------------------------------------#
# Now for FastCMH

    def test_create_flags_fastcmh1(self):
        '''Wrapper has flags for fastcmh:
           use_intervals=True, use_covariate=True
        '''
        sig = createSigPatSearch(use_intervals=True, use_covariate=True)
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchFastCmh"
        self.assertEqual(sigclass, solution)


    def test_create_flags_fastcmh2(self):
        '''Wrapper has flags for fastcmh:
           use_combinations=False, use_covariate=True
        '''
        sig = createSigPatSearch(use_combinations=False, use_covariate=True)
        sigclass = type(sig).__name__
        solution = "_SignificantIntervalSearchFastCmh"
        self.assertEqual(sigclass, solution)


#----------------------------------------------------------------------#
# Now for FACS


    def test_create_flags_facs1(self):
        '''Wrapper has flags for facs:
           use_intervals=False, use_covariate=False (default)
        '''
        sig = createSigPatSearch(use_intervals=False)
        sigclass = type(sig).__name__
        solution = "_SignificantItemsetSearchFacs"
        self.assertEqual(sigclass, solution)


    def test_create_flags_facs2(self):
        '''Wrapper has flags for facs:
           use_intervals=False, use_covariate=False
        '''
        sig = createSigPatSearch(use_intervals=False, use_covariate=False)
        sigclass = type(sig).__name__
        solution = "_SignificantItemsetSearchFacs"
        self.assertEqual(sigclass, solution)


    def test_create_flags_facs3(self):
        '''Wrapper has flags for facs:
           use_intervals=False, use_covariate=True
        '''
        sig = createSigPatSearch(use_intervals=False, use_covariate=True)
        sigclass = type(sig).__name__
        solution = "_SignificantItemsetSearchFacs"
        self.assertEqual(sigclass, solution)


    def test_create_flags_facs4(self):
        '''Wrapper has flags for facs:
           use_combinations=True, use_covariate=False (default)
        '''
        sig = createSigPatSearch(use_combinations=True)
        sigclass = type(sig).__name__
        solution = "_SignificantItemsetSearchFacs"
        self.assertEqual(sigclass, solution)


    def test_create_flags_facs5(self):
        '''Wrapper has flags for facs:
           use_combinations=True, use_covariate=False 
        '''
        sig = createSigPatSearch(use_combinations=True, use_covariate=False)
        sigclass = type(sig).__name__
        solution = "_SignificantItemsetSearchFacs"
        self.assertEqual(sigclass, solution)


    def test_create_flags_facs6(self):
        '''Wrapper has flags for facs:
           use_combinations=True, use_covariate=True
        '''
        sig = createSigPatSearch(use_combinations=True, use_covariate=True)
        sigclass = type(sig).__name__
        solution = "_SignificantItemsetSearchFacs"
        self.assertEqual(sigclass, solution)


#----------------------------------------------------------------------#
#force error messages with bad flags

    def test_errors_from_bad_flags1(self):
        '''No flags
           Throws an error, no parameters set
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch()


    def test_errors_from_bad_flags2(self):
        '''No flags for intervals, combinations or method,
           only for covariate (which is not enough)
           Throws an error, no parameters set
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch(use_covariate=True)


    def test_errors_from_bad_flags3(self):
        '''When one boolan flag is not a boolean
           Throws an error
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch(use_combinations="x", use_covariate=True)


    def test_errors_from_bad_flags4(self):
        '''When use_intervals is None
           Throws an error
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch(use_intervals=None, use_covariate=True)


    def test_errors_from_bad_flags5(self):
        '''When use_combinations is None
           Throws an error
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch(use_combinations=None, use_covariate=True)


    def test_errors_from_bad_flags6(self):
        '''When use_combinations is not boolean
           Throws an error
        '''
        with self.assertRaises(ValueError):
            sig = createSigPatSearch(use_combinations=2.3, use_covariate=True)

#----------------------------------------------------------------------#
# check alpha and lmax values

    def test_alpha_and_lmax1(self):
        '''Check default alpha value
        '''
        sig = createSigPatSearch(use_combinations=True)
        default_alpha = 0.05
        self.assertEqual(sig.get_alpha(), default_alpha)


    def test_alpha_and_lmax2(self):
        '''Check default lmax value
        '''
        sig = createSigPatSearch(use_combinations=True)
        default_lmax = 0
        self.assertEqual(sig.get_lmax(), default_lmax)


    def test_alpha_and_lmax3(self):
        '''Setting alpha value
        '''
        sig = createSigPatSearch(use_combinations=True)
        new_alpha = 0.05
        sig.set_alpha(new_alpha)
        self.assertEqual(sig.get_alpha(), new_alpha)


    def test_alpha_and_lmax4(self):
        '''Setting lmax value
        '''
        sig = createSigPatSearch(use_combinations=True)
        new_lmax = 10
        sig.set_lmax(new_lmax)
        self.assertEqual(sig.get_lmax(), new_lmax)



# end of tests
#----------------------------------------------------------------------#





if __name__ == '__main__':
    unittest.main()












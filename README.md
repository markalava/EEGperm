# EEGperm

This repository contains R code used in the article:

> Wheldon, M. C., Anderson, M. J., and Johnson, B. W. (2007) "Identifying Treatment Effects in Multi-Channel Measurements in Electroencephalographic Studies: Multivariate Permutation Tests and Multiple Comparisons", *Australia and New Zealand Journal of Statistics,* 49(4) 397-413.  https://doi.org/10.1111/j.1467-842X.2007.00491.x

*R* code and examples for the permutation *t*- and *F*-tests are in the 'R/' directory in individual files. Each file contains the source code to run the tests followed by a toy example with randomly generated data. You should be able to source each script as-is to get example re-creations of Figures 1 and 2 from the manuscript. The objects generated by the functions can be examined for the test statistics, quantiles, and adjusted p-values.

Please see the file 'LICENSE' in this directory for the licence and terms of use. 

Please cite the original paper and this repository if you make use of these methods and/or the code. Note, though, that the functions `binary.v` and `perm.t.test3` in 'R/t_tests_Section_3.1_Functions_and_Example.R' (the permuation t-tests) are taken from Venables and Ripley (2002, Sect 5.7), slightly modified in the latter case, so recommend checking and citing that reference if you use those functions. 

#### References
Venables and Ripley (2002) *Modern Applied Statistics with S*, 4th edn. New York: Springer

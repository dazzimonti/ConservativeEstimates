# ConservativeEstimates 0.2.0

## Major changes

* Added a `NEWS.md` file to track changes to the package.

* Online choice of the number of active dimensions: the number of active dimensions can now be chosen online. Given a range of possible values for _q_, the function `selectQdims` selects at the same time _q_ and the active dimensions.

* The order of the inputs in `ProbaMax` and `ProbaMin` is changed. They also can take as input a range of possible values for _q_ instead of a single value _q_.

* The function `conservativeEstimate` now uses a range of possible _q_ to initialize the probability estimates.

## Bug fixes

* The functions `ANMC_Gauss` and `MC_Gauss` are now more efficient as we do not keep all inner simulations in memory at the same time. 

* Minor changes in `ANMC_Gauss`: mStar is computed with an equivalent more efficient formula that prevents possible memory overflows. Some for loops in the initial parameter estimation phase were fine tuned for efficiency and to prevent possible memory overflows.




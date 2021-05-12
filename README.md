Fit linear models in Futhark with a robustness towards ill-conditioned problems.
Model parameters are estimated by solving the linear least squares
equations based on a pivoting QR-decomposition.
It is well suited for small, numerically instable problems.

The QR-decomposition is _rank revealing_, meaning we obtain the rank
of the regressor matrix, `X`. The rank is used to ensure a meaningful
fit: if the problem is ill-conditioned so that (at least for numerical
work) the regressor matrix is rank deficient, parameters will
be dropped accordingly until this is no longer the case.
It should be possible to extend this to also output the condition number of `X`,
or at least an estimate of this---untested sample code can be found in `lm.fut`.

This package seeks to match the output of `lm` in the R language.
Therefore, the QR-decomposition and parts of the least squares solver
is a 1:1 translation of sequential FORTRAN routines from the famous
LINPACK library (see `linpack_d`).

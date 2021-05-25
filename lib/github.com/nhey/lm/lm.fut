-- | Linear models using rank-revealing QR decomposition.
-- Linear dependencies and near singularities are dealt
-- with by dropping parameters.
-- This module is desirable when the model is fit on a small
-- number of observations relative to the number of parameters
-- and high accuracy is needed.

module lm_f64 = {
  import "../../diku-dk/linalg/linalg"
  import "linpack"

  module T = f64
  type real = T.t

  module linalg = mk_linalg T
  module linpack = linpack_d

  -- Dotproduct ignoring nans.
  let dotprod_nan [n] (xs: [n]real) (ys: [n]real): real =
    reduce (+) 0 (map2 (\x y -> if T.isnan y then 0 else x*y) xs ys)

  -- Back substitution for solving upper triangular linear systems.
  -- Floating point nan values in `U` are considered padding.
  let back_substitution [n] (U: [n][n]real) (y: [n]real): [n]real =
    let x = replicate n (T.i64 0)
    in loop x for j in 0..<n do
      let i = n - j - 1
      let sumx = dotprod_nan x[i+1:n] U[i,i+1:n]
      let x[i] = (y[i] T.- sumx) T./ U[i,i]
      in x

  let identity (n: i64): [n][n]real =
    tabulate_2d n n (\i j ->if j == i then T.i64 1 else T.i64 0)

  -- Given an upper triangular matrix `U`, compute `(U^T U)^{-1}`.
  -- Also returns the intermediate result `U^{-1}`.
  -- * If fed transpose of `L` from Cholesky decomposition of a symmetric,
  -- positive definite square matrix `A`, result is `A^{-1}`.
  -- * If fed `R` from QR decomposition of `X`, result is `(X^T X)^{-1}`
  -- since `X^T X = R^T Q^T Q R = R^T R`.
  -- The name is an homage to similar functionality in the R language.
  let chol2inv [n] (U: [n][n]real): ([n][n]real, [n][n]real) =
    let UinvT = map (back_substitution U) (identity n)
    let Uinv = transpose UinvT
    -- Compute `(U^T U)^{-1} = U^{-1} (U^T)^{-1} = U^{-1} (U^{-1})^T`.
    in (linalg.matmul Uinv UinvT, Uinv)

  -- Estimate the condition number of `X` given `R` from
  -- QR decomposition of `X` as `k(X) = k(R) = ||R|| ||R^{-1}||`.
  -- The estimate is correct within a factor 2, roughly.
  -- TODO is this not exact when using rank-revealing QR?
  -- TODO with rank-revealing QR, compute as max abs diag / min abs diag?
  -- TODO expose this somehow
  let cond_est [n] (R: [n][n]real) (Rinv: [n][n]real): real =
    let frobenius_norm_sq A = flatten A |> \x -> linalg.dotprod x x
    in T.sqrt <| (frobenius_norm_sq R) T.* (frobenius_norm_sq Rinv)

  -- Extract lower triangular matrix by setting values above
  -- main diagonal to NAN. Main diagonal is determined by
  -- the rank.
  -- If fed output from `dqrdc2`, it will extract transposed `R`
  -- from `QR`, seing as `dqrdc2` output itself is transposed.
  let lower_triangular_nan [m][n] (rank: i64) (L: [m][n]real): [m][n]real =
    map2 (\j ->
            map2 (\i ele -> if i+1 > rank || j+1 > rank
                            then T.nan
                            else ele
                 ) (iota n)
         ) (iota m) L

  type results [p] = { params: [p]real, cov_params: [p][p]real, rank: i64 }

  let fit [n][p] (X': [p][n]real) (y: [n]real): results [p] =
    let (A', pivot, qraux, rank) = linpack.dqrdc2 (copy X') 1e-7
    -- The shared dimension of `Q` and `R` is `k = min(n,p,r)`
    -- where `r` is the rank of `X`. Fitting a linear regression
    -- with `p` parameters requires at least `p` datapoints, so
    -- we always have `n >= p`. Moreover the rank of `X` is at
    -- most `p`. Thus it is always the case that `k = r`.
    let R = transpose (lower_triangular_nan rank A'[:p,:p]) -- k x k
    let (cov_params, Rinv) = chol2inv R
    -- Find least squares solution to `Xb = y`. Substituting
    -- `X = QR` into the LLS equations `X^T X b = X^T y`, we get
    -- `((QR)^T QR) b = (QR)^T y <=> (R^T Q^T QR) b = R^T Q^T y`.
    -- Now since `Q^T Q = I`, we have `R^T R b = R^T Q^T y`.
    -- "Since R is full rank", we can ignore the `R^T` factor
    -- yielding `R b = Q^T y`.
    -- With `R` being upper triangular, this last equation may be
    -- solved using back substitution. But we already have
    -- `R^{-1}` and premultiplying this is less sequential.
    let Q'y = linpack.dqrqty A' qraux rank y
    let Q'y = Q'y[:p]
    let beta = linalg.matvecmul_row Rinv Q'y
    -- Pivot to match original `X`.
    let beta = scatter (replicate p (T.i64 0)) pivot beta
    let pp = p*p
    let pivot_2d = map (\i -> map (\j -> (i,j)) pivot) pivot
                   |> flatten :> [pp](i64,i64)
    let cov_params = scatter_2d (replicate pp (T.i64 0) |> unflatten p p)
                                pivot_2d (cov_params |> flatten :> [pp]real)
    -- Set padding values equal to zero.
    let cov_params = map (map (\x -> if T.isnan x then 0 else x)) cov_params
    in { params = beta, cov_params = cov_params, rank = rank }
}

module linpack_d = {
  module T = f64
  type t = f64

  -- Compute `Q'y` given output from `dqrdc2` (`Q` is reconstructed
  -- from `x` and `qraux`).
  let dqrqty [n][p] (x: [p][n]t) (qraux: [p]t) (k: i64) (y: [n]t) =
    let ju = i64.min k (n-1)
    in loop (qty) = (copy y) for j < ju do
      if qraux[j] T.== (T.i64 0)
      then qty
      else let t = - qraux[j] * qty[j]
           let t =
             loop (t) for i0 < (n - j - 1) do -- for i in [j+1..<n] do
               let i = i0 + j + 1
               in t - x[j,i] * qty[i]
           let t = t / qraux[j]
           let qty[j] = qty[j] + t * qraux[j]
           let qty =
             loop (qty) for i0 < (n - j - 1) do
               let i = i0 + j + 1
               in qty with [i] = qty[i] + t * x[j,i]
           in qty


  local let dotprod xs ys : t =
    T.(reduce (+) (i64 0) (map2 (*) xs ys))

  local let dnrm2 xs =
    T.sqrt (dotprod xs xs)

  -- Set the sign of scalar a to match the sign of scalar b.
  -- Undefined when b is zero.
  local let dsign a b =
    T.((sgn b) * (abs a))

  let dqrdc2 [n][p] (x: *[p][n]t) (tol: t): ([p][n]t, [p]i64, [p]t, i64) =
    let qraux = replicate p (T.i64 0)
    let work = replicate (2*p) (T.i64 0) |> unflatten 2 p
    let jpvt = iota p
    -- compute the norms of the columns of x.
    let (qraux, work) =
      loop (qraux, work) for j < p do
        let nrm = dnrm2 x[j,:]
        let qraux[j] = nrm
        let work[0,j] = nrm
        let work[1,j] = if T.(nrm == i64 0) then T.i64 1 else nrm
        in (qraux, work)
    -- perform the householder reduction of x.
    let lup = i64.min n p
    let k = p + 1
    --   cycle the columns from l to p left-to-right until one
    --   with non-negligible norm is located.  a column is considered
    --   to have become negligible if its norm has fallen below
    --   tol times its original norm.  the check for l .le. k
    --   avoids infinite cycling.
    let (x, jpvt, qraux, _, k) =
      loop (x, jpvt, qraux, work, k) for l < lup do
        let (x, jpvt, qraux, work, k) =
          loop (x, jpvt, qraux, work, k) = (x, jpvt, qraux, work, k)
            while (l+1 < k) && T.(qraux[l] < work[1,l] * tol) do
              let lp1 = l + 1
              let x =
                loop (x) for i < n do
                  let t = x[l,i]
                  let x =
                    loop (x) for j0 < (p - lp1) do
                      let j = j0 + lp1
                      let x[j-1,i] = x[j,i]
                      in  x
                  let x[p-1, i] = t
                  in  x
              let i = jpvt[l]
              let t   = qraux[l]
              let tt  = work[0,l]
              let ttt = work[1,l]

              let (jpvt,qraux,work) =
                loop (jpvt,qraux,work)
                for j0 < (p - lp1) do
                  let j = j0 + lp1
                  let jpvt[j-1] = jpvt[j]
                  let qraux[j-1] = qraux[j]
                  let work[0,j-1] = work[0,j]
                  let work[1,j-1] = work[1,j]
                  in  (jpvt,qraux,work)

              let jpvt[p-1] = i
              let qraux[p-1] = t
              let work[0,p-1] = tt
              let work[1,p-1] = ttt
              let k = k - 1
              in  (x, jpvt, qraux, work, k)
        -- compute the householder transformation for column l.
        let (qraux, work, x) =
          if l+1 == n
          then (qraux, work, x)
          else let nrmxl = dnrm2 x[l,l:]
               in if T.(nrmxl == i64 0)
                  then (qraux, work, x)
                  else let nrmxl = if T.(x[l,l] != i64 0)
                                   then dsign nrmxl x[l,l]
                                   else nrmxl
                       let x =
                         loop (x) for i0 < n - l do
                           let i = i0 + l
                           in x with [l,i] = x[l,i] / nrmxl
                       let x[l,l] = T.(i64 1 + x[l,l])
                       -- apply the transformation to the remaining columns,
                       -- updating the norms.
                       let (qraux, work, x) =
                         loop (qraux, work, x) for j in (l+1..<p) do
                           let t =
                             loop (t) = (T.i64 0) for i0 < n - l do
                               let i = i0 + l
                               in t - x[l,i] * x[j,i]
                           let t = t / x[l,l]
                           let x =
                             loop (x) for i0 < n - l do
                               let i = i0 + l
                               in x with [j,i] = x[j,i] + t * x[l,i]
                           in if T.(qraux[j] == i64 0)
                              then (qraux, work, x)
                              else let tt = T.(((abs x[j,l])/qraux[j]) ** i64 2)
                                   let tt = T.(i64 1 - tt)
                                   let tt = T.(max tt (i64 0))
                                   let t = tt
                                   in if T.(abs t >= f64 1e-6)
                                      then let qraux[j] = T.(qraux[j] * sqrt t)
                                           in (qraux, work, x)
                                      else let qraux[j] = dnrm2 x[j,l+1:]
                                           let work[0,j] = qraux[j]
                                           in (qraux, work, x)
                       let qraux[l] = x[l,l]
                       let x[l,l] = T.(i64 (-1) * nrmxl)
                       in (qraux, work, x)
        in (x, jpvt, qraux, work, k)
    let k = i64.min (k-1) n
    in (x, jpvt, qraux, k)
}

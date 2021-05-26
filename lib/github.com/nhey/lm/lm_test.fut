-- | ignore

import "lm"

module lm = lm_f64

-- Test various fits. First data set exhibits full rank.
-- Second data set is rank deficient, but does not trigger pivoting.
-- Third data set is rank deficient and will require pivoting.
-- ==
-- entry: fit
-- input { [[1.0f64, 86.0f64], [1.0f64, 76.0f64], [1.0f64, 92.0f64],
--          [1.0f64, 90.0f64], [1.0f64, 86.0f64], [1.0f64, 84.0f64],
--          [1.0f64, 93.0f64], [1.0f64, 100.0f64], [1.0f64, 87.0f64],
--          [1.0f64, 86.0f64], [1.0f64, 74.0f64], [1.0f64, 98.0f64],
--          [1.0f64, 97.0f64], [1.0f64, 84.0f64], [1.0f64, 91.0f64],
--          [1.0f64, 34.0f64], [1.0f64, 45.0f64], [1.0f64, 56.0f64],
--          [1.0f64, 44.0f64], [1.0f64, 82.0f64], [1.0f64, 72.0f64],
--          [1.0f64, 55.0f64], [1.0f64, 71.0f64], [1.0f64, 50.0f64],
--          [1.0f64, 23.0f64], [1.0f64, 39.0f64], [1.0f64, 28.0f64],
--          [1.0f64, 32.0f64], [1.0f64, 22.0f64], [1.0f64, 25.0f64],
--          [1.0f64, 29.0f64], [1.0f64, 7.0f64],  [1.0f64, 26.0f64],
--          [1.0f64, 19.0f64], [1.0f64, 15.0f64], [1.0f64, 20.0f64],
--          [1.0f64, 26.0f64], [1.0f64, 28.0f64], [1.0f64, 17.0f64],
--          [1.0f64, 22.0f64], [1.0f64, 30.0f64], [1.0f64, 25.0f64],
--          [1.0f64, 20.0f64], [1.0f64, 47.0f64], [1.0f64, 32.0f64]]
--         [62.0f64, 72.0f64, 75.0f64, 55.0f64, 64.0f64, 21.0f64,
--          64.0f64, 80.0f64, 67.0f64, 72.0f64, 42.0f64, 76.0f64,
--          76.0f64, 41.0f64, 48.0f64, 76.0f64, 53.0f64, 60.0f64,
--          42.0f64, 78.0f64, 29.0f64, 48.0f64, 55.0f64, 29.0f64,
--          21.0f64, 47.0f64, 81.0f64, 36.0f64, 22.0f64, 44.0f64,
--          15.0f64,  7.0f64, 42.0f64,  9.0f64, 21.0f64, 21.0f64,
--          16.0f64, 16.0f64,  9.0f64, 14.0f64, 12.0f64, 17.0f64,
--           7.0f64, 34.0f64,  8.0f64] }
-- output { [10.60349832f64, 0.59485944f64]
--          [[9.30974511e-02f64, -1.34857729e-03f64],
--           [-1.34857729e-03f64, 2.56600331e-05f64]]
--          2i64 }
-- input { [[1.0f64, 1.0f64,   0.5f64,       0.866025f64,  0.866025f64,  0.5f64,
--            1.0f64,  0.0f64],
--          [1.0f64, 2.0f64,   0.866025f64,  0.5f64,       0.866025f64, -0.5f64,
--            0.0f64, -1.0f64],
--          [1.0f64, 4.0f64,   0.866025f64, -0.5f64,      -0.866025f64, -0.5f64,
--           -0.0f64,  1.0f64],
--          [1.0f64, 5.0f64,   0.5f64,      -0.866025f64, -0.866025f64,  0.5f64,
--            1.0f64,  0.0f64],
--          [1.0f64, 7.0f64,  -0.5f64,      -0.866025f64,  0.866025f64,  0.5f64,
--           -1.0f64, -0.0f64],
--          [1.0f64, 9.0f64,  -1.0f64,      -0.0f64,       0.0f64,      -1.0f64,
--            1.0f64,  0.0f64],
--          [1.0f64, 14.0f64,  0.866025f64,  0.5f64,       0.866025f64, -0.5f64,
--            0.0f64, -1.0f64],
--          [1.0f64, 17.0f64,  0.5f64,      -0.866025f64, -0.866025f64,  0.5f64,
--            1.0f64,  0.0f64]]
--         [4074.389219f64, 6842.004064f64, 6378.106146f64, 5561.734381f64,
--          6480.266660f64, 4951.876764f64, 4043.239055f64, 7192.433631f64] }
-- output { [ 4894.79236538f64,  -48.66940662f64,  988.41905285f64,
--           -2798.41799501f64, 1187.96283884f64, -903.18463376f64,
--             580.34347733f64,    0.0f64]
--          [[ 2.29786725e+00f64,  2.40563102e-02f64, -1.51490260e+00f64,
--             3.79778089e+00f64, -3.03688247e+00f64,  1.62654830e+00f64,
--            -2.24802784e+00, 0.0f64],
--           [ 2.40563102e-02f64,  6.94444444e-03f64, -4.30288848e-02f64,
--             1.32808172e-01f64, -9.27143845e-02f64,  5.69177737e-02f64,
--            -7.26674213e-02, 0.0f64],
--           [-1.51490260e+00f64, -4.30288848e-02f64,  1.36242131e+00f64,
--            -2.93294738e+00f64,  2.26439826e+00f64, -1.27130628e+00f64,
--             1.63601192e+00, 0.0f64],
--           [ 3.79778089e+00f64,  1.32808172e-01f64, -2.93294738e+00f64,
--             8.13796208e+00f64, -6.18817882e+00f64,  3.42623488e+00f64,
--            -4.49976695e+00, 0.0f64],
--           [-3.03688247e+00f64, -9.27143845e-02f64,  2.26439826e+00f64,
--            -6.18817882e+00f64,  4.96987773e+00f64, -2.58127309e+00f64,
--             3.55443710e+00, 0.0f64],
--           [ 1.62654830e+00f64,  5.69177737e-02f64, -1.27130628e+00f64,
--             3.42623488e+00f64, -2.58127309e+00f64,  1.81908015e+00f64,
--            -1.90043540e+00, 0.0f64],
--           [-2.24802784e+00f64, -7.26674213e-02f64,  1.63601192e+00f64,
--            -4.49976695e+00f64,  3.55443710e+00f64, -1.90043540e+00f64,
--             2.81624398e+00, 0.0f64],
--           [ 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64,
--             0.0f64]]
--          7i64 }
-- input { [[1929f64,  2804f64,    0f64,  4836f64,  3595f64],
--          [2804f64,  6386f64,    0f64, 10174f64,  6737f64],
--          [   0f64,     0f64,    0f64,     0f64,     0f64],
--          [4836f64, 10174f64,    0f64, 21229f64, 14415f64],
--          [3595f64,  6737f64,    0f64, 14415f64, 10420f64]]
--         [6.369616873e-01f64, 2.697867138e-01f64, 4.097352394e-02f64,
--          1.652763553e-02f64, 8.132702392e-01f64] }
-- output { [ 0.00011335f64, 0.00028335f64, 0.0f64,
--           -0.00103702f64, 0.00129035f64]
--          [[ 8.06407657e-06f64, -4.34739617e-06f64,  0.00000000e+00f64,
--             5.79730793e-06f64, 0.00000000e+00f64],
--           [-4.34739617e-06f64,  2.87097032e-06f64,  0.00000000e+00f64,
--            -3.72025866e-06f64, 0.00000000e+00f64],
--           [ 0.00000000e+00f64,  0.00000000e+00f64,  0.00000000e+00f64,
--             0.00000000e+00f64, 0.00000000e+00f64],
--           [ 5.79730793e-06f64, -3.72025866e-06f64,  0.00000000e+00f64,
--             5.96059073e-06f64, 0.00000000e+00f64],
--           [ 0.00000000e+00f64,  0.00000000e+00f64,  0.00000000e+00f64,
--             0.00000000e+00f64, 0.00000000e+00f64]]
--          4i64 }
entry fit [n][p] (X: [n][p]lm.real) (y: [n]lm.real) =
  let res = lm.fit (transpose X) y
  in (res.params, res.cov_params, res.rank)

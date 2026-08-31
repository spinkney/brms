  /*
   * The semi-orthogonal Householder-reflector transform below is adapted
   * from Seth Axen's stan_semiorthogonal_transforms:
   * https://github.com/sethaxen/stan_semiorthogonal_transforms
   * upstream revision 8f47b5760a896e363ee5e15966c4f96e017c42eb.
   *
   * Modifications for brms namespace the functions, separate the Gaussian
   * prior from the pure transform, validate dimensions, correct the final
   * special-orthogonal reflector, and add the centered Helmert lift.
   *
   * MIT License
   *
   * Copyright (c) 2022 Seth Axen <seth@sethaxen.com> and contributors
   *
   * Permission is hereby granted, free of charge, to any person obtaining a
   * copy of this software and associated documentation files (the "Software"),
   * to deal in the Software without restriction, including without limitation
   * the rights to use, copy, modify, merge, publish, distribute, sublicense,
   * and/or sell copies of the Software, and to permit persons to whom the
   * Software is furnished to do so, subject to the following conditions:
   *
   * The above copyright notice and this permission notice shall be included
   * in all copies or substantial portions of the Software.
   *
   * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   * DEALINGS IN THE SOFTWARE.
   */

  // Apply the elementary Householder reflector I - tau * v * v'.
  matrix fm_apply_reflector_brms(real tau, vector v, matrix B) {
    if (tau == 0) {
      return B;
    }
    return B - (tau * v) * (v' * B);
  }

  // The determinant of a product of Householder reflectors.
  real fm_determinant_from_tau_brms(vector tau) {
    int num_negative_det = 0;
    for (j in 1:num_elements(tau)) {
      if (tau[j] != 0) {
        num_negative_det += 1;
      }
    }
    if (fmod(num_negative_det, 2) == 0) {
      return 1;
    }
    return -1;
  }

  // Stable positive-diagonal variant of LAPACK's DLARFGP.
  tuple(real, vector, real) fm_get_reflector_brms(real x1, vector xtail) {
    real tau;
    vector[num_elements(xtail)] v;
    real beta;
    real xtail_norm = norm2(xtail);
    if (xtail_norm == 0) {
      tau = 2 * (x1 <= 0);
      v = xtail;
      beta = abs(x1);
      return (tau, v, beta);
    }
    {
      real xnorm = hypot(x1, xtail_norm);
      real eta;
      beta = -xnorm;
      if (x1 < 0) {
        beta = -beta;
      }
      if (beta >= 0) {
        eta = x1 - beta;
      } else {
        real gamma;
        beta = -beta;
        gamma = x1 + beta;
        eta = -xtail_norm * (xtail_norm / gamma);
      }
      tau = -eta / beta;
      v = xtail / eta;
    }
    return (tau, v, beta);
  }

  tuple(real, real) fm_get_reflector_brms(real x1) {
    real tau = 2 * (x1 <= 0);
    real beta = abs(x1);
    return (tau, beta);
  }

  int fm_semiorthogonal_num_params_brms(int N, int K, int special) {
    if (N < 1 || K < 1 || K > N) {
      reject("A semi-orthogonal matrix must have 1 <= K <= N.");
    }
    if (special != 0 && special != 1) {
      reject("Argument 'special' must be either zero or one.");
    }
    return N * K - (K * (K - 1)) %/% 2
           - (N == K && special == 1);
  }

  // Convert Gaussian reflector coordinates to LAPACK-style (V, tau) factors.
  tuple(matrix, vector) fm_reflector_factors_brms(
      vector y, int N, int K, int special) {
    matrix[N, K] factors = rep_matrix(0, N, K);
    vector[K] tau = zeros_vector(K);
    int Nlower = N - 1;
    int iy = 1;
    int maxcols = K - (special == 1 && N == K);
    int expected = fm_semiorthogonal_num_params_brms(N, K, special);
    if (num_elements(y) != expected) {
      reject("Wrong number of semi-orthogonal reflector coordinates; found ",
             num_elements(y), ", expected ", expected, ".");
    }
    if (maxcols > 0) {
      for (j in 1:maxcols) {
        real x1 = y[iy];
        if (Nlower > 0) {
          vector[Nlower] xtail = segment(y, iy + 1, Nlower);
          tuple(real, vector[Nlower], real) reflector =
            fm_get_reflector_brms(x1, xtail);
          tau[j] = reflector.1;
          factors[(j + 1):N, j] = reflector.2;
        } else {
          tuple(real, real) reflector = fm_get_reflector_brms(x1);
          tau[j] = reflector.1;
        }
        factors[j, j] = 1;
        iy += Nlower + 1;
        Nlower -= 1;
      }
    }
    if (special == 1 && N == K) {
      factors[N, N] = 1;
      if (N == 1) {
        tau[N] = 0;
      } else {
        // The scalar final reflector has determinant 1 - tau[N].
        tau[N] = 1 - fm_determinant_from_tau_brms(head(tau, N - 1));
      }
    }
    return (factors, tau);
  }

  // Form the thin Q matrix from LAPACK-style Householder factors.
  matrix fm_factors_to_Q_brms(matrix factors, vector tau) {
    int N = rows(factors);
    int K = cols(factors);
    matrix[N, K] Q = rep_matrix(0, N, K);
    int Nlower = N - K + 1;
    for (j in reverse(linspaced_int_array(K, 1, K))) {
      vector[Nlower] v = factors[j:N, j];
      Q[j, j] = 1;
      Q[j:N, j:K] = fm_apply_reflector_brms(
        tau[j], v, Q[j:N, j:K]
      );
      Nlower += 1;
    }
    return Q;
  }

  matrix fm_semiorthogonal_constrain_brms(
      vector y, int N, int K, int special) {
    tuple(matrix[N, K], vector[K]) Q_fact =
      fm_reflector_factors_brms(y, N, K, special);
    return fm_factors_to_Q_brms(Q_fact.1, Q_fact.2);
  }

  // Isometric Helmert lift from R^(L-1) into the zero-sum subspace of R^L.
  vector fm_sum_to_zero_constrain_brms(vector y) {
    int N = num_elements(y);
    vector[N + 1] z = zeros_vector(N + 1);
    real sum_w = 0;
    for (ii in 1:N) {
      int i = N - ii + 1;
      real w = y[i] * inv_sqrt(i * (i + 1.0));
      sum_w += w;
      z[i] += sum_w;
      z[i + 1] -= i * w;
    }
    return z;
  }

  // A semi-orthogonal frame whose columns also sum exactly to zero.
  matrix fm_centered_semiorthogonal_constrain_brms(
      vector y, int L, int K, int special) {
    matrix[L - 1, K] Q_small = fm_semiorthogonal_constrain_brms(
      y, L - 1, K, special
    );
    matrix[L, K] Q;
    for (k in 1:K) {
      Q[, k] = fm_sum_to_zero_constrain_brms(Q_small[, k]);
    }
    return Q;
  }

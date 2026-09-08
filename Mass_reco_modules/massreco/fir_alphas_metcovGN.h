// Gauss-Newton solver for multiplicative mass constraint with MET covariance.
// Assumes the nonlinearity is on system A: M_A = mvisA * sqrt((1+alpha1)*(1+alpha2))
// and system B mass is fixed (mvisB). If you need the reverse, swap A/B inputs.
//
// Uses same tiny linear algebra helpers as before.

#include <array>
#include <cmath>
#include <iostream>
#include "3nu.h"

inline vec2 add(const vec2 &a,const vec2 &b){ return vec2{a[0]+b[0], a[1]+b[1]}; }
inline vec2 subv(const vec2 &a,const vec2 &b){ return vec2{a[0]-b[0], a[1]-b[1]}; }
inline vec2 scalev(double s,const vec2 &v){ return vec2{s*v[0], s*v[1]}; }

// Solve H x = b for 2x2 H using inv2
inline bool solve2(const mat2 &H, const vec2 &b, vec2 &x) {
    bool ok;
    mat2 Hi = inv2(H, ok);
    if (!ok) return false;
    x = mat2_vec(Hi, b);
    return true;
}

/*
 multiplicative Gauss-Newton solver.

 Inputs:
  k1,k2 : 2D vectors (columns of U)
  met   : measured MET (2D)
  Sigma_met : MET covariance (2x2)
  mvisA  : visible mass for system A (the multiplicative side)
  mvisB  : visible mass for system B (fixed)
  sigmaM_A, sigmaM_B : mass uncertainties (sigma, not variance)

 Params:
  maxIter (default 20), tol for delta norm

 Outputs:
  alpha_star (2), cov_alpha_star (2x2), chi2_mass (value of lambda * f^2/(1+... denom) approx),
  returns true if converged / success, false otherwise.

 Note: this routine does NOT enforce alpha >= 0. If you need non-negativity, enforce after solution (or modify iteration).
*/
bool FitAlphasMultiplicativeGN(
    TVector3 &Lxy1, TVector3 &Lxy2,
    TVector3 &MET,
    const mat2 &Sigma_met,
    double mvisA, double mvisB,
    vec2 &alpha_star,
    mat2 &cov_alpha_star,
    double &chi2_mass,
    vec2 &alpha,
    int maxIter = 20,
    double tol = 1e-6
) {
    // Build U
    mat2 U;
    U[0][0] = Lxy1.X(); U[0][1] = Lxy2.X();
    U[1][0] = Lxy1.Y(); U[1][1] = Lxy2.Y();
    
    vec2 met;
	met[0] = MET.X(); met[1] = MET.Y();
	
    bool ok;
    mat2 Uinv = inv2(U, ok);
    if (!ok) return false; // nearly collinear
    // alpha_hat (linear solution)
    vec2 alpha_hat = mat2_vec(Uinv, met);

    // Sigma_alpha and its inverse W
    mat2 tmp = mat2_mul(Uinv, Sigma_met);
    mat2 UinvT = mat2_transpose(Uinv);
    mat2 Sigma_alpha = mat2_mul(tmp, UinvT);

    mat2 W = inv2(Sigma_alpha, ok);
    if (!ok) return false;
	double sigmaM_A = 0.02*mvisA, sigmaM_B = 0.02*mvisB;
    // lambda weight for mass constraint
    double lambda = 1.0 / (sigmaM_A*sigmaM_A + sigmaM_B*sigmaM_B);

    // Initialize alpha to alpha_hat as a starting guess
    //vec2 alpha = alpha_hat;
    
    // Iterative Gauss-Newton
    for (int iter = 0; iter < maxIter; ++iter) {
        // compute M_A(alpha) = mvisA * sqrt((1+alpha1)*(1+alpha2))
        double onep_a1 = 1.0 + alpha[0];
        double onep_a2 = 1.0 + alpha[1];
        if (onep_a1 <= 0 || onep_a2 <= 0) {
            // invalid region -> cannot take sqrt; abort
            //cout<<"BAD SQRT"<<endl;
            return false;
        }
        double MA = mvisA * std::sqrt(onep_a1 * onep_a2);
        double MB = mvisB; // fixed in this case

        double f = MA - MB; // scalar residual

        // Jacobian J = df/dalpha (1x2)
        // dMA/dalpha1 = (1/2) * mvisA * sqrt((1+alpha2)/(1+alpha1)) = MA / (2*(1+alpha1))
        // similarly for alpha2
        double dMA_da1 = MA / (2.0 * onep_a1);
        double dMA_da2 = MA / (2.0 * onep_a2);
        vec2 J = { dMA_da1, dMA_da2 }; // row vector

        // Build gradient: g = W*(alpha - alpha_hat) + lambda * J^T * f
        vec2 alpha_diff = subv(alpha, alpha_hat);
        vec2 Wad = mat2_vec(W, alpha_diff); // W*(alpha - alpha_hat)
        vec2 JT_f = scalev(f, J);          // J^T * f (componentwise)
        vec2 grad = add(Wad, scalev(lambda, JT_f)); // 2-vector

        // Build approximate Hessian: H = W + lambda * J^T J  (2x2)
        mat2 JTJ = outer(J, J); // 2x2
        mat2 H = mat2_mul(W, mat2{{{1,0},{0,1}}}); // copy W
        for (int a=0;a<2;++a) for (int b=0;b<2;++b)
            H[a][b] += lambda * JTJ[a][b];

        // Solve H * delta = -grad
        vec2 rhs = scalev(-1.0, grad);
        vec2 delta;
        if (!solve2(H, rhs, delta)) {
            // ill-conditioned Hessian
			cout<<"BAD HESSIAN"<<endl;
            return false;
        }

        // update
        alpha = add(alpha, delta);

        // convergence check
        double dn = std::sqrt(delta[0]*delta[0] + delta[1]*delta[1]);
        if (dn < tol) break;
    }

    // final outputs
    alpha_star = alpha;

    // Final covariance approximation = H^{-1}
    // Recompute final J and H to return covariance
    double onep_a1f = 1.0 + alpha_star[0];
    double onep_a2f = 1.0 + alpha_star[1];
    if (onep_a1f <= 0 || onep_a2f <= 0) return false;
    double MAf = mvisA * std::sqrt(onep_a1f * onep_a2f);
    vec2 Jf = { MAf / (2.0 * onep_a1f), MAf / (2.0 * onep_a2f) };
    mat2 JTJf = outer(Jf, Jf);
    mat2 Hfinal = W;
    for (int a=0;a<2;++a) for (int b=0;b<2;++b)
        Hfinal[a][b] += lambda * JTJf[a][b];

    mat2 Hinv = inv2(Hfinal, ok);
    if (!ok) return false;
    cov_alpha_star = Hinv;

    // chi2 contribution from mass term: approx 0.5 * lambda * f^2 ; but for reporting follow same as earlier:
    // We'll compute chi2_mass = 0.5 * lambda * f^2 (or lambda * f^2 / 2)
    double MA_final = MAf;
    double MB_final = mvisB;
    double ffinal = MA_final - MB_final;
    //chi2_mass = 0.5 * lambda * ffinal * ffinal;
	    
    double beta_eff = dot(Jf, mat2_vec(Sigma_alpha, Jf));
    chi2_mass = lambda * ffinal*ffinal / (1 + lambda*beta_eff);
    
    return true;
}

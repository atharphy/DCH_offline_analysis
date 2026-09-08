// fit_alphas_metcov.h  (or paste into a .cpp)
// Implements the closed-form α* solution with MET covariance (see user's PDF). :contentReference[oaicite:1]{index=1}
#pragma once
#include <array>
#include <cmath>
#include <iostream>


using vec2 = std::array<double,2>;
using mat2 = std::array<std::array<double,2>,2>;

// --- small linear-algebra helpers for 2x2 ---
inline double det2(const mat2 &A) {
    return A[0][0]*A[1][1] - A[0][1]*A[1][0];
}
inline mat2 inv2(const mat2 &A, bool &ok) {
    double d = det2(A);
    const double eps = 1e-12;
    ok = std::fabs(d) > eps;
    if (!ok) return mat2{{{0,0},{0,0}}};
    double id = 1.0/d;
    mat2 R;
    R[0][0] =  A[1][1]*id;
    R[0][1] = -A[0][1]*id;
    R[1][0] = -A[1][0]*id;
    R[1][1] =  A[0][0]*id;
    return R;
}
inline vec2 mat2_vec(const mat2 &A, const vec2 &v) {
    return vec2{ A[0][0]*v[0] + A[0][1]*v[1],
                 A[1][0]*v[0] + A[1][1]*v[1] };
}
inline mat2 mat2_mul(const mat2 &A, const mat2 &B) {
    mat2 C{};
    for (int i=0;i<2;++i) for (int j=0;j<2;++j) {
        C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j];
    }
    return C;
}
inline mat2 mat2_transpose(const mat2 &A){
    mat2 B;
    B[0][0]=A[0][0]; B[0][1]=A[1][0];
    B[1][0]=A[0][1]; B[1][1]=A[1][1];
    return B;
}
inline vec2 vec2_add(const vec2 &a, const vec2 &b){ return vec2{a[0]+b[0], a[1]+b[1]}; }
inline vec2 vec2_sub(const vec2 &a, const vec2 &b){ return vec2{a[0]-b[0], a[1]-b[1]}; }
inline vec2 scalar_vec(double s, const vec2 &v){ return vec2{ s*v[0], s*v[1] }; }
inline mat2 outer(const vec2 &a, const vec2 &b){
    mat2 M; M[0][0]=a[0]*b[0]; M[0][1]=a[0]*b[1]; M[1][0]=a[1]*b[0]; M[1][1]=a[1]*b[1]; return M;
}
inline double dot(const vec2 &a, const vec2 &b){ return a[0]*b[0] + a[1]*b[1]; }

// Add small diagonal regularization (if desired) to a matrix:
inline mat2 add_diag(const mat2 &A, double eps) {
    mat2 B = A;
    B[0][0] += eps;
    B[1][1] += eps;
    return B;
}

inline double quad(TVector3 &a, const mat2 &M, TVector3 &b) {
    return a.X()*(M[0][0]*b.X() + M[0][1]*b.Y()) +
           a.Y()*(M[1][0]*b.X() + M[1][1]*b.Y());
}

/*
 Solve for alpha in MET = alpha * L using MET covariance.

 Inputs:
   L         : lepton transverse vector (px, py)
   MET       : measured MET vector
   Sigma_met : 2x2 MET covariance

 Outputs:
   alpha     : best-fit alpha
   sigma2    : variance of alpha
   chi2min   : minimum chi^2

 Returns:
   true on success, false if covariance is singular or L too small
*/
bool SolveOneNeutrinoAlpha(
    TVector3 &L,
    TVector3 &MET,
    const mat2 &Sigma_met,
    double &alpha,
    double &sigma2,
    double &chi2min
) {
    // Check lepton magnitude
    /*double L2 = L.X()+L.X() + L.Y()+L.Y();
    if (L2 < 1e-12) return false;

    // Invert MET covariance
    bool ok;
    mat2 W = inv2(Sigma_met, ok);
    if (!ok) return false;

    // Compute numerator and denominator
    double num = quad(L, W, MET); // L^T W MET
    double den = quad(L, W, L);   // L^T W L

    if (den < 1e-14) return false;

    // Best-fit alpha
    alpha = num / den;

    // Variance
    sigma2 = 1.0 / den;

    // Minimum chi2
    double METW_MET = quad(MET, W, MET);
    chi2min = METW_MET - (num*num)/den;*/

    // Invert MET covariance
    bool ok;
    mat2 W = inv2(Sigma_met, ok);

    // Best-fit alpha
    alpha = MET.Mag();//L.Mag();

    // Variance
    sigma2 = 1.0;
    // Minimum chi2
    chi2min = 0;
    
    return true;
}

/*
 Returns:
   - alpha_star (size 2)
   - cov_alpha_star (2x2)
   - chi2min (minimized chi2)
   - M2star (profiled M^2)
 Returns true on success, false on failure (matrix singular).
 Notes:
   - m2_ll_A, m2_ll_B are the dilepton mass squared (GeV^2 or whatever units you use).
   - sigmaM_A, sigmaM_B are the uncertainties (σ_M) for the masses (same units as sqrt(m2) or here used as sigma in eq; the formula uses sigma^2 so we square them).
   - No positivity constraint on alpha is enforced here.
*/
bool FitAlphasWithMETCov(
    TVector3 &Lxy1, TVector3 &Lxy2,
    TVector3 &MET,
    const mat2 &Sigma_met,
    double m2_ll_A, double m2_ll_B,
    double alpha_A, double alpha_B,
    vec2 &alpha_star,
	mat2 &cov_alpha_star,
	double &chi2min, 
	double &M2star
) {
	double Delta = 0000.0; // GeV^2 (tune this!)

    // Build U = [k1 k2] as 2x2 (columns are k1,k2)
    mat2 U;
    U[0][0] = Lxy1.X(); U[0][1] = Lxy2.X();
    U[1][0] = Lxy1.Y(); U[1][1] = Lxy2.Y();
	
	vec2 alpha = {{alpha_A, alpha_B}};
	
	vec2 met;
	met[0] = MET.X(); met[1] = MET.Y();
	
    // Invert U
    bool ok;
    mat2 Uinv = inv2(U, ok);
    if (!ok) {
        // nearly collinear k1/k2 -> cannot invert U
        return false;
    }

    // Step 1: alpha_hat = U^{-1} * met
    vec2 alpha_hat = mat2_vec(Uinv, met);

    // Sigma_alpha = Uinv * Sigma_met * (Uinv)^T
    mat2 tmp = mat2_mul(Uinv, Sigma_met);
    mat2 UinvT = mat2_transpose(Uinv);
    mat2 Sigma_alpha = mat2_mul(tmp, UinvT);

    // W = Sigma_alpha^{-1}
    mat2 W = inv2(Sigma_alpha, ok);
    if (!ok) {
        // extremely ill-conditioned Sigma_alpha
        return false;
    }
    
	double sigmaM_A = 0.02*sqrt(m2_ll_A), sigmaM_B = 0.02*sqrt(m2_ll_B);
    // Step 2: mass-constraint width lambda = 1 / (sigmaM_A^2 + sigmaM_B^2)
    double lambda = 1.0 / (sigmaM_A*sigmaM_A + sigmaM_B*sigmaM_B);

    // c vector: [ m2_ll_A, - m2_ll_B ]  (so c^T alpha + d = m2_A*alpha1 - m2_B*alpha2 + d)
    vec2 c = { m2_ll_A, -m2_ll_B };
    double d = m2_ll_A - m2_ll_B;

    // beta = c^T Sigma_alpha c  (scalar)
    // compute Sigma_alpha * c
    vec2 Sigma_alpha_c = mat2_vec(Sigma_alpha, c);
    double beta = dot(c, Sigma_alpha_c);

    // compute scalar s = c^T alpha_hat + d
    double s = dot(c, alpha_hat) + d;

    // denominator denom = 1 + lambda * beta
    double denom = 1.0 + lambda * beta;
    if (std::fabs(denom) < 1e-15) {
        // numerical trouble
        return false;
    }

   /* // Step 3: alpha* = alpha_hat - (lambda * Sigma_alpha * c * s) / denom
    vec2 term = scalar_vec(lambda * s / denom, Sigma_alpha_c); // Sigma_alpha*c scaled
    
    alpha_star = vec2{ alpha_hat[0] - term[0], alpha_hat[1] - term[1] };

    // Sigma_alpha,* = Sigma_alpha - (lambda * Sigma_alpha * c c^T Sigma_alpha) / denom
    // Implement numerator matrix N = lambda * (Sigma_alpha * c) (Sigma_alpha * c)^T
    mat2 numerator = outer(Sigma_alpha_c, Sigma_alpha_c);
    mat2 correction;
    double scale = lambda / denom;
    for (int i=0;i<2;++i) for (int j=0;j<2;++j)
        correction[i][j] = scale * numerator[i][j];
        
    for (int i=0;i<2;++i) for (int j=0;j<2;++j)
        cov_alpha_star[i][j] = Sigma_alpha[i][j] - correction[i][j];

    // Step 4: minimized chi2: chi2min = lambda * s^2 / denom
    chi2min = 0.00001*lambda * s * s / denom;*/
	
	double s_raw = dot(c, alpha_hat) + d;
	double s_eff = 0.0;
	if (std::abs(s_raw) > Delta)
    	s_eff = (s_raw > 0 ? 1.0 : -1.0) * (std::abs(s_raw) - Delta);
	// Step 3: alpha* = alpha_hat - (lambda * Sigma_alpha * c * s) / denom
    vec2 term = scalar_vec(lambda * s_eff / denom, Sigma_alpha_c);
    
    alpha_star = vec2{ alpha_hat[0] - term[0], alpha_hat[1] - term[1] };

    // Sigma_alpha,* = Sigma_alpha - (lambda * Sigma_alpha * c c^T Sigma_alpha) / denom
    // Implement numerator matrix N = lambda * (Sigma_alpha * c) (Sigma_alpha * c)^T
    mat2 numerator = outer(Sigma_alpha_c, Sigma_alpha_c);
    mat2 correction;
    double scale = lambda / denom;
    for (int i=0;i<2;++i) for (int j=0;j<2;++j)
        correction[i][j] = scale * numerator[i][j];
        
    for (int i=0;i<2;++i) for (int j=0;j<2;++j)
        cov_alpha_star[i][j] = Sigma_alpha[i][j] - correction[i][j];

    // Step 4: minimized chi2: chi2min = lambda * s^2 / denom
   chi2min = 0.00001*lambda * s_eff * s_eff / denom;



    // Optional: profiled common mass squared M2* per eq (6)
    // XA = m2_ll_A * (1 + alpha1*)
    // XB = m2_ll_B * (1 + alpha2*)
    
    double XA = m2_ll_A * (1.0 + alpha_star[0]);
    double XB = m2_ll_B * (1.0 + alpha_star[1]);
    double wA = 1.0 / (sigmaM_A*sigmaM_A);
    double wB = 1.0 / (sigmaM_B*sigmaM_B);
    M2star = (XA * wA + XB * wB) / (wA + wB);
	
    return true;
}

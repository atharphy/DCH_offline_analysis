#pragma once
#include <array>
#include <cmath>
#include <iostream>

using vec2 = std::array<double,2>;
using mat2 = std::array<std::array<double,2>,2>;

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

bool SolveOneNeutrinoAlpha(
    TVector3 &L,
    TVector3 &MET,
    const mat2 &Sigma_met,
    double &alpha,
    double &sigma2,
    double &chi2min
) {

    bool ok;
    mat2 W = inv2(Sigma_met, ok);

    alpha = MET.Mag();

    sigma2 = 1.0;

    chi2min = 0;

    return true;
}

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
	double Delta = 0000.0;

    mat2 U;
    U[0][0] = Lxy1.X(); U[0][1] = Lxy2.X();
    U[1][0] = Lxy1.Y(); U[1][1] = Lxy2.Y();

	vec2 alpha = {{alpha_A, alpha_B}};

	vec2 met;
	met[0] = MET.X(); met[1] = MET.Y();

    bool ok;
    mat2 Uinv = inv2(U, ok);
    if (!ok) {

        return false;
    }

    vec2 alpha_hat = mat2_vec(Uinv, met);

    mat2 tmp = mat2_mul(Uinv, Sigma_met);
    mat2 UinvT = mat2_transpose(Uinv);
    mat2 Sigma_alpha = mat2_mul(tmp, UinvT);

    mat2 W = inv2(Sigma_alpha, ok);
    if (!ok) {

        return false;
    }

	double sigmaM_A = 0.02*sqrt(m2_ll_A), sigmaM_B = 0.02*sqrt(m2_ll_B);

    double lambda = 1.0 / (sigmaM_A*sigmaM_A + sigmaM_B*sigmaM_B);

    vec2 c = { m2_ll_A, -m2_ll_B };
    double d = m2_ll_A - m2_ll_B;

    vec2 Sigma_alpha_c = mat2_vec(Sigma_alpha, c);
    double beta = dot(c, Sigma_alpha_c);

    double s = dot(c, alpha_hat) + d;

    double denom = 1.0 + lambda * beta;
    if (std::fabs(denom) < 1e-15) {

        return false;
    }

	double s_raw = dot(c, alpha_hat) + d;
	double s_eff = 0.0;
	if (std::abs(s_raw) > Delta)
    	s_eff = (s_raw > 0 ? 1.0 : -1.0) * (std::abs(s_raw) - Delta);

    vec2 term = scalar_vec(lambda * s_eff / denom, Sigma_alpha_c);

    alpha_star = vec2{ alpha_hat[0] - term[0], alpha_hat[1] - term[1] };

    mat2 numerator = outer(Sigma_alpha_c, Sigma_alpha_c);
    mat2 correction;
    double scale = lambda / denom;
    for (int i=0;i<2;++i) for (int j=0;j<2;++j)
        correction[i][j] = scale * numerator[i][j];

    for (int i=0;i<2;++i) for (int j=0;j<2;++j)
        cov_alpha_star[i][j] = Sigma_alpha[i][j] - correction[i][j];

   chi2min = 0.00001*lambda * s_eff * s_eff / denom;

    double XA = m2_ll_A * (1.0 + alpha_star[0]);
    double XB = m2_ll_B * (1.0 + alpha_star[1]);
    double wA = 1.0 / (sigmaM_A*sigmaM_A);
    double wB = 1.0 / (sigmaM_B*sigmaM_B);
    M2star = (XA * wA + XB * wB) / (wA + wB);

    return true;
}

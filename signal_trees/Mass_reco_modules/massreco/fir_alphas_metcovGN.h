#include <array>
#include <cmath>
#include <iostream>
#include "3nu.h"

inline vec2 add(const vec2 &a,const vec2 &b){ return vec2{a[0]+b[0], a[1]+b[1]}; }
inline vec2 subv(const vec2 &a,const vec2 &b){ return vec2{a[0]-b[0], a[1]-b[1]}; }
inline vec2 scalev(double s,const vec2 &v){ return vec2{s*v[0], s*v[1]}; }

inline bool solve2(const mat2 &H, const vec2 &b, vec2 &x) {
    bool ok;
    mat2 Hi = inv2(H, ok);
    if (!ok) return false;
    x = mat2_vec(Hi, b);
    return true;
}

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

    mat2 U;
    U[0][0] = Lxy1.X(); U[0][1] = Lxy2.X();
    U[1][0] = Lxy1.Y(); U[1][1] = Lxy2.Y();

    vec2 met;
	met[0] = MET.X(); met[1] = MET.Y();

    bool ok;
    mat2 Uinv = inv2(U, ok);
    if (!ok) return false;

    vec2 alpha_hat = mat2_vec(Uinv, met);

    mat2 tmp = mat2_mul(Uinv, Sigma_met);
    mat2 UinvT = mat2_transpose(Uinv);
    mat2 Sigma_alpha = mat2_mul(tmp, UinvT);

    mat2 W = inv2(Sigma_alpha, ok);
    if (!ok) return false;
	double sigmaM_A = 0.02*mvisA, sigmaM_B = 0.02*mvisB;

    double lambda = 1.0 / (sigmaM_A*sigmaM_A + sigmaM_B*sigmaM_B);

    for (int iter = 0; iter < maxIter; ++iter) {

        double onep_a1 = 1.0 + alpha[0];
        double onep_a2 = 1.0 + alpha[1];
        if (onep_a1 <= 0 || onep_a2 <= 0) {

            return false;
        }
        double MA = mvisA * std::sqrt(onep_a1 * onep_a2);
        double MB = mvisB;

        double f = MA - MB;

        double dMA_da1 = MA / (2.0 * onep_a1);
        double dMA_da2 = MA / (2.0 * onep_a2);
        vec2 J = { dMA_da1, dMA_da2 };

        vec2 alpha_diff = subv(alpha, alpha_hat);
        vec2 Wad = mat2_vec(W, alpha_diff);
        vec2 JT_f = scalev(f, J);
        vec2 grad = add(Wad, scalev(lambda, JT_f));

        mat2 JTJ = outer(J, J);
        mat2 H = mat2_mul(W, mat2{{{1,0},{0,1}}});
        for (int a=0;a<2;++a) for (int b=0;b<2;++b)
            H[a][b] += lambda * JTJ[a][b];

        vec2 rhs = scalev(-1.0, grad);
        vec2 delta;
        if (!solve2(H, rhs, delta)) {

			cout<<"BAD HESSIAN"<<endl;
            return false;
        }

        alpha = add(alpha, delta);

        double dn = std::sqrt(delta[0]*delta[0] + delta[1]*delta[1]);
        if (dn < tol) break;
    }

    alpha_star = alpha;

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

    double MA_final = MAf;
    double MB_final = mvisB;
    double ffinal = MA_final - MB_final;

    double beta_eff = dot(Jf, mat2_vec(Sigma_alpha, Jf));
    chi2_mass = lambda * ffinal*ffinal / (1 + lambda*beta_eff);

    return true;
}

/*
bessel_ref.h

テスト用の整数次ベッセル関数 J_l(x) / 変形ベッセル関数 K_l(x)。

C++17 の std::cyl_bessel_j / std::cyl_bessel_k は libc++ (macOS/AppleClang) が
未実装のため、全プラットフォームで同一結果になる自前実装を使う
(LP モード分散方程式の求解に必要な l <= 2, 0 < x <= ~10 の範囲を対象)。

- J_l : 冪級数 (この範囲では倍精度で完全収束)
- K_0/K_1 : Abramowitz & Stegun 9.8.5-9.8.8 の多項式近似 (|誤差| < 1e-7)、
  高次は上方漸化式 K_{n+1} = K_{n-1} + (2n/x) K_n (K は上方向に安定)
- I_0/I_1 : A&S 9.8.1/9.8.3 (K の小引数枝でのみ使用、x <= 2 < 3.75)

精度検証: libstdc++ の std::cyl_bessel_j/k と l=0..2, x=0.05..10 で比較し
相対誤差 J < 1e-14 / K < 2e-7 (tests 追加時のコミットに記録)。
分散方程式の二分法には十分な精度。
*/
#ifndef TESTS_BESSEL_REF_H_
#define TESTS_BESSEL_REF_H_

#include <cmath>

// J_l(x) : 第 1 種ベッセル関数 (冪級数)。x <= ~10, l <= 3 を想定
static inline double bessel_ref_j(int l, double x)
{
    double term = 1.0;                       // (x/2)^l / l!
    for (int k = 1; k <= l; k++) term *= (x / 2.0) / k;
    double sum = term;
    const double q = x * x / 4.0;
    for (int m = 1; m < 80; m++) {
        term *= -q / (double(m) * double(m + l));
        sum += term;
        if (std::fabs(term) < 1e-18 * std::fabs(sum)) break;
    }
    return sum;
}

// I_0(x) : A&S 9.8.1 (|x| <= 3.75)
static inline double bessel_ref_i0_small(double x)
{
    const double t = (x / 3.75) * (x / 3.75);
    return 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492
              + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
}

// I_1(x) : A&S 9.8.3 (|x| <= 3.75)
static inline double bessel_ref_i1_small(double x)
{
    const double t = (x / 3.75) * (x / 3.75);
    return x * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934
             + t * (0.02658733 + t * (0.00301532 + t * 0.00032411))))));
}

// K_l(x) : 第 2 種変形ベッセル関数 (x > 0)
static inline double bessel_ref_k(int l, double x)
{
    double k0, k1;
    if (x <= 2.0) {
        // A&S 9.8.5 / 9.8.7
        const double y = x * x / 4.0;
        k0 = -std::log(x / 2.0) * bessel_ref_i0_small(x)
           + (-0.57721566 + y * (0.42278420 + y * (0.23069756
              + y * (0.03488590 + y * (0.00262698 + y * (0.00010750
              + y * 0.00000740))))));
        k1 = std::log(x / 2.0) * bessel_ref_i1_small(x)
           + (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579
              + y * (-0.18156897 + y * (-0.01919402 + y * (-0.00110404
              + y * (-0.00004686)))))));
    } else {
        // A&S 9.8.6 / 9.8.8
        const double y = 2.0 / x;
        const double c = std::exp(-x) / std::sqrt(x);
        k0 = c * (1.25331414 + y * (-0.07832358 + y * (0.02189568
           + y * (-0.01062446 + y * (0.00587872 + y * (-0.00251540
           + y * 0.00053208))))));
        k1 = c * (1.25331414 + y * (0.23498619 + y * (-0.03655620
           + y * (0.01504268 + y * (-0.00780353 + y * (0.00325614
           + y * (-0.00068245)))))));
    }
    if (l == 0) return k0;
    // 上方漸化式 (K は上方向に安定)
    double km1 = k0, kn = k1;
    for (int n = 1; n < l; n++) {
        const double kp1 = km1 + (2.0 * n / x) * kn;
        km1 = kn;
        kn = kp1;
    }
    return kn;
}

#endif  // TESTS_BESSEL_REF_H_

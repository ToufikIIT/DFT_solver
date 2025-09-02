#include <bits/stdc++.h>
using namespace std;

const int N_grid = 2000; 
const double r_min = 1e-6;
const double r_max = 50.0;
const double dr = (r_max - r_min) / (N_grid - 1);
const double Z = 2.0; 
const int ell = 0;   

vector<double> setup_r_grid()
{
    vector<double> r(N_grid);
    for (int i = 0; i < N_grid; ++i)
        r[i] = r_min + i * dr;
    return r;
}

double trapz(const vector<double> &f)
{
    double s = 0.0;
    for (int i = 1; i < (int)f.size(); ++i)
        s += 0.5 * (f[i - 1] + f[i]) * dr;
    return s;
}

double trapz_radial(const vector<double> &f, const vector<double> &r)
{
    double s = 0.0;
    for (int i = 1; i < (int)f.size(); ++i)
    {
        double w1 = 4.0 * M_PI * r[i - 1] * r[i - 1];
        double w2 = 4.0 * M_PI * r[i] * r[i];
        s += 0.5 * (f[i - 1] * w1 + f[i] * w2) * dr;
    }
    return s;
}

double trapz_weighted_radial(const vector<double> &f, const vector<double> &g, const vector<double> &r)
{
    double s = 0.0;
    for (int i = 1; i < (int)f.size(); ++i)
    {
        double w1 = 4.0 * M_PI * r[i - 1] * r[i - 1];
        double w2 = 4.0 * M_PI * r[i] * r[i];
        s += 0.5 * (f[i - 1] * g[i - 1] * w1 + f[i] * g[i] * w2) * dr;
    }
    return s;
}

vector<double> calc_v_x(const vector<double> &n)
{
    vector<double> vx(N_grid, 0.0);
    const double c = pow(3.0 / M_PI, 1.0 / 3.0);
    const double small = 1e-20;
    for (int i = 0; i < N_grid; ++i)
    {
        double ni = n[i];
        if (ni <= small)
            vx[i] = 0.0;
        else
            vx[i] = -c * pow(ni, 1.0 / 3.0);
    }
    return vx;
}

double exchange_energy(const vector<double> &n, const vector<double> &r)
{
    const double C = -0.75 * pow(3.0 / M_PI, 1.0 / 3.0);
    vector<double> e(n.size(), 0.0);
    const double small = 1e-20;
    for (int i = 0; i < N_grid; ++i)
    {
        double ni = n[i];
        if (ni <= small)
            e[i] = 0.0;
        else
            e[i] = C * pow(ni, 4.0 / 3.0);
    }
    return trapz_radial(e, r);
}

vector<double> calc_hartree(const vector<double> &n, const vector<double> &r)
{
    vector<double> VH(N_grid, 0.0);
    vector<double> A(N_grid, 0.0), B(N_grid, 0.0);

    vector<double> n_r2(N_grid), n_r(N_grid);
    for (int i = 0; i < N_grid; ++i)
    {
        n_r2[i] = n[i] * r[i] * r[i];
        n_r[i] = n[i] * r[i];
    }

    A[0] = 0.0;
    for (int i = 1; i < N_grid; ++i)
    {
        A[i] = A[i - 1] + 2.0 * M_PI * (n_r2[i - 1] + n_r2[i]) * dr;
    }
    B[N_grid - 1] = 0.0;
    for (int i = N_grid - 2; i >= 0; --i)
    {
        B[i] = B[i + 1] + 2.0 * M_PI * (n_r[i + 1] + n_r[i]) * dr;
    }

    for (int i = 0; i < N_grid; ++i)
    {
        VH[i] = A[i] / r[i] + B[i];
    }
    return VH;
}
vector<double> build_v_eff(const vector<double> &r,
                           const vector<double> &VH,
                           const vector<double> &vxc)
{
    vector<double> V(N_grid, 0.0);
    for (int i = 0; i < N_grid; ++i)
    {
        // centrifugal barrier is handled separately in k(r) inside Numerov
        V[i] = -Z / r[i] + VH[i] + vxc[i];
    }
    return V;
}
void numerov_solve(const vector<double> &r, const vector<double> &Veff, double E, vector<double> &u)
{
    const int N = r.size();
    u.assign(N, 0.0);

    const double ell_d = static_cast<double>(ell);
    auto k = [&](int i) -> double
    {
        double centrifugal = (ell_d * (ell_d + 1.0)) / (r[i] * r[i]);
        return 2.0 * (E - Veff[i]) - centrifugal;
    };
    u[0] = 0.0;
    u[1] = (r[1] - r[0]); 
    const double h2 = dr * dr;
    for (int i = 1; i < N - 1; ++i)
    {
        double k_im1 = k(i - 1);
        double k_i = k(i);
        double k_ip1 = k(i + 1);

        double c_ip1 = 1.0 + h2 * k_ip1 / 12.0;

        u[i + 1] = ((2.0 * (1.0 - 5.0 * h2 * k_i / 12.0)) * u[i] - (1.0 + h2 * k_im1 / 12.0) * u[i - 1]) / c_ip1;
        if (!isfinite(u[i + 1]))
            u[i + 1] = 0.0;
    }
}

void normalize_u(vector<double> &u)
{
    vector<double> u2(u.size());
    for (int i = 0; i < u.size(); ++i)
        u2[i] = u[i] * u[i];
    double nrm = sqrt(trapz(u2));
    if (nrm > 0.0)
        for (auto &x : u)
            x /= nrm;
}

vector<double> density_from_u(const vector<double> &r, const vector<double> &u, double occ = 2.0)
{
    vector<double> n(N_grid, 0.0);
    const double inv4pi = 1.0 / (4.0 * M_PI);
    for (int i = 0; i < N_grid; ++i)
    {
        double ri = r[i];
        double denom = max(ri * ri, 1e-20);
        n[i] = occ * (u[i] * u[i]) * inv4pi / denom;
    }
    return n;
}

double find_ground_energy(const vector<double> &r, const vector<double> &Veff,
                          double Emin = -3.0, double Emax = -0.1, int itmax = 60)
{
    vector<double> uL, uR, uM;
    numerov_solve(r, Veff, Emin, uL);
    numerov_solve(r, Veff, Emax, uR);
    double fL = uL.back();
    double fR = uR.back();
    if (fL * fR > 0.0)
    {
        double expand = 1.5;
        int tries = 0;
        while (fL * fR > 0.0 && tries < 10)
        {
            Emin *= expand;
            Emax *= (1.0 / expand);
            numerov_solve(r, Veff, Emin, uL);
            numerov_solve(r, Veff, Emax, uR);
            fL = uL.back();
            fR = uR.back();
            tries++;
        }
    }

    double a = Emin, b = Emax;
    for (int it = 0; it < itmax; ++it)
    {
        double m = 0.5 * (a + b);
        numerov_solve(r, Veff, m, uM);
        double fm = uM.back();
        if (fL * fm <= 0.0)
        {
            b = m;
            fR = fm;
        }
        else
        {
            a = m;
            fL = fm;
        }
        if (fabs(b - a) < 1e-9)
            break;
    }
    return 0.5 * (a + b);
}

int main()
{
    auto r = setup_r_grid();

    vector<double> u(N_grid, 0.0);
    for (int i = 0; i < N_grid; ++i)
    {
        // R_1s(r) ~ 2 Z^{3/2} e^{-Z r}; u=r R
        double Ri = 2.0 * pow(Z, 1.5) * exp(-Z * r[i]);
        u[i] = r[i] * Ri;
    }
    normalize_u(u);

    vector<double> n = density_from_u(r, u, 2.0);

    const int max_iter = 200;
    const double mix_alpha = 0.5; 
    const double E_tol = 1e-7;
    const double dens_tol = 1e-6;

    double Etot_prev = 1e100;

    for (int iter = 1; iter <= max_iter; ++iter)
    {
        // Hartree + XC
        vector<double> VH = calc_hartree(n, r);
        vector<double> vxc = calc_v_x(n);

        // Effective potential
        vector<double> Veff = build_v_eff(r, VH, vxc);

        // Orbital energy by shooting
        double eps = find_ground_energy(r, Veff);

        numerov_solve(r, Veff, eps, u);
        normalize_u(u);

        // New density
        vector<double> n_new = density_from_u(r, u, 2.0);

        double max_rel = 0.0;
        for (int i = 0; i < N_grid; ++i)
        {
            double mixed = mix_alpha * n_new[i] + (1.0 - mix_alpha) * n[i];
            max_rel = max(max_rel, fabs(mixed - n[i]) / (1e-8 + n[i]));
            n[i] = mixed;
        }

        // Rebuild Hartree/XC with mixed density for energy evaluation
        VH = calc_hartree(n, r);
        vxc = calc_v_x(n);

        // Energies
        // KS sum of eigenvalues (two electrons in same orbital)
        double Esum = 2.0 * eps;

        // Hartree double-counting correction: -1/2 ∫ n VH d^3r
        double EH_dc = -0.5 * trapz_weighted_radial(n, VH, r);

        // Exchange energy and potential contribution: E_x - ∫ n v_x d^3r
        double Ex = exchange_energy(n, r);
        double int_nvx = trapz_weighted_radial(n, vxc, r);

        double Etot = Esum + EH_dc + Ex - int_nvx;

        cout << "Iter " << setw(3) << iter
             << "  eps = " << setw(12) << setprecision(8) << eps
             << "  Etot = " << setw(12) << setprecision(8) << Etot
             << "  densΔ(max rel) = " << scientific << max_rel << "\n";

        if (fabs(Etot - Etot_prev) < E_tol && max_rel < dens_tol)
        {
            cout << "\nConverged:\n";
            cout << "  Orbital energy (1s): " << setprecision(10) << eps << " Ha\n";
            cout << "  Total KS energy    : " << setprecision(10) << Etot << " Ha\n";
            break;
        }
        Etot_prev = Etot;

        if (iter == max_iter)
        {
            cout << "\nReached max iterations.\n";
            cout << "  Orbital energy (1s): " << setprecision(10) << eps << " Ha\n";
            cout << "  Total KS energy    : " << setprecision(10) << Etot << " Ha\n";
        }
    }

    return 0;
}

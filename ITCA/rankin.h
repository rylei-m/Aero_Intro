#include <iostream>
#include <iomanip>
#include <cmath>
#include <optional>

using namespace std;

struct RankineResult {
    double x, psi, u, v, cp;
};

optional<RankineResult> rankine_oval(double Q, double y) {
    double x0 = -Q + sqrt(Q*Q + 1.0);
    double denom = tan(y / Q);
    if (fabs(denom) < 1e-12) return nullopt;

    double x2 = x0 + x0 - y*y + 2.0*y*x0/denom;
    if (x2 < 0.0) return nullopt;

    double x = sqrt(x2);
    double psi = y + Q*atan2(y, x + x0) - Q*atan2(y, x - x0);

    double rplus  = (x + x0)*(x + x0) + y*y;
    double rminus = (x - x0)*(x - x0) + y*y;
    double u = 1.0 + Q*(x + x0)/rplus - Q*(x - x0)/rminus;
    double v = Q*y/rplus - Q*y/rminus;
    double cp = 1.0 - u*u - v*v;

    return RankineResult{x, psi, u, v, cp};
}

int main() {
    cout << "\n=== PROGRAM RANKIN (Rankine Oval) ===\n";
    double Q;
    cout << "DIMENSIONLESS SOURCE STRENGTH = ";
    if (!(cin >> Q)) return 0;

    cout << setw(10) << "X"
         << setw(10) << "PSI"
         << setw(10) << "U"
         << setw(10) << "V"
         << setw(10) << "CP\n";

    while (true) {
        cout << "Y/A = ";
        double y; if (!(cin >> y)) break;
        auto res = rankine_oval(Q, y);
        if (!res) {
            cout << "Y is too large (no real solution).\n";
            continue;
        }
        cout << setw(10) << fixed << setprecision(3) << res->x
             << setw(10) << res->psi
             << setw(10) << res->u
             << setw(10) << res->v
             << setw(10) << res->cp << "\n";
    }
    return 0;
}

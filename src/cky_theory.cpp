// Canonical estimators in CKY (2026), Definitions: CATE Estimators,
// Tree Construction, and Sample Splitting and Estimators.
// No package defaults, RNG, pruning, variance penalties, or ancestor fallback.
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
struct Stats {
    int n = 0, n1 = 0;
    long double s0 = 0, s1 = 0, z = 0;
    void add(double y, int d, double xi) {
        ++n; n1 += d;
        if (d) { s1 += y; z += static_cast<long double>(y) / xi; }
        else { s0 += y; z -= static_cast<long double>(y) / (1 - xi); }
    }
};
Stats difference(const Stats& total, const Stats& left) {
    Stats r;
    r.n = total.n - left.n; r.n1 = total.n1 - left.n1;
    r.s0 = total.s0 - left.s0; r.s1 = total.s1 - left.s1;
    r.z = total.z - left.z;
    return r;
}
long double dim_mean(const Stats& s) {
    return s.s1 / s.n1 - s.s0 / (s.n - s.n1);
}
double estimate(const Stats& s, int rule) {
    if (rule == 0) return s.n ? static_cast<double>(s.z / s.n) : 0.0;
    if (s.n1 == 0 || s.n1 == s.n) return 0.0;
    return static_cast<double>(dim_mean(s));
}
struct Node {
    int left = -1, right = -1, feature = -1, depth = 0;
    double cut = NA_REAL, gain = NA_REAL, value = 0;
    Stats train, est;
};
struct Builder {
    const double *x, *y;
    const int *d;
    int n, p, maxdepth, rule;
    double xi;
    std::vector<Node> nodes;
    int grow(const std::vector<int>& rows, int depth) {
        const int id = static_cast<int>(nodes.size());
        nodes.emplace_back();
        nodes[id].depth = depth;
        Stats total;
        for (int i : rows) total.add(y[i], d[i], xi);
        nodes[id].train = total;
        if (depth >= maxdepth || rows.size() < 2) return id;
        if (rule != 0 && (total.n1 < 2 || total.n - total.n1 < 2)) return id;

        long double best = -1;
        int feature = -1;
        double cut = NA_REAL;
        for (int j = 0; j < p; ++j) {
            std::vector<int> order(rows);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return x[a + n*j] < x[b + n*j] ||
                    (x[a + n*j] == x[b + n*j] && a < b);
            });
            Stats left;
            for (size_t k = 0; k + 1 < order.size(); ++k) {
                const int row = order[k];
                left.add(y[row], d[row], xi);
                if (x[row + n*j] == x[order[k+1] + n*j]) continue;
                const Stats right = difference(total, left);
                if (rule != 0 && (left.n1 == 0 || left.n1 == left.n ||
                                  right.n1 == 0 || right.n1 == right.n)) continue;
                long double score;
                if (rule == 0 || rule == 1) {
                    const long double delta = rule == 0 ?
                        left.z / left.n - right.z / right.n :
                        dim_mean(left) - dim_mean(right);
                    score = static_cast<long double>(left.n) * right.n / total.n * delta * delta;
                } else {
                    // Exactly the profiled OLS SSE gain: separate arm means.
                    const long double delta1 = left.s1 / left.n1 - right.s1 / right.n1;
                    const int l0 = left.n - left.n1, r0 = right.n - right.n1;
                    const long double delta0 = left.s0 / l0 - right.s0 / r0;
                    score = static_cast<long double>(left.n1) * right.n1 / total.n1 * delta1 * delta1 +
                        static_cast<long double>(l0) * r0 / (total.n - total.n1) * delta0 * delta0;
                }
                if (!std::isfinite(score)) throw std::runtime_error("Nonfinite split objective");
                // Deterministic tie order: first coordinate, then smallest cutoff.
                // A zero maximum remains a valid maximizer; no gain-based pruning.
                const long double tolerance = 64 * std::numeric_limits<double>::epsilon() *
                    std::max(std::abs(best), std::abs(score));
                if (best < 0 || score > best + tolerance) {
                    best = score; feature = j; cut = x[row + n*j];
                }
            }
        }
        if (feature < 0) return id;
        std::vector<int> leftrows, rightrows;
        for (int i : rows) (x[i + n*feature] <= cut ? leftrows : rightrows).push_back(i);
        nodes[id].feature = feature; nodes[id].cut = cut;
        nodes[id].gain = static_cast<double>(best);
        const int left = grow(leftrows, depth + 1);
        const int right = grow(rightrows, depth + 1);
        nodes[id].left = left; nodes[id].right = right;
        return id;
    }
};
}

extern "C" SEXP cky_fit(SEXP tx, SEXP ty, SEXP td, SEXP ex, SEXP ey, SEXP ed,
                         SEXP rule_, SEXP xi_, SEXP depth_) {
    try {
        const int n = Rf_nrows(tx), p = Rf_ncols(tx), ne = Rf_nrows(ex);
        Builder builder{REAL(tx), REAL(ty), INTEGER(td), n, p,
                        Rf_asInteger(depth_), Rf_asInteger(rule_), Rf_asReal(xi_), {}};
        std::vector<int> rows(n);
        std::iota(rows.begin(), rows.end(), 0);
        builder.grow(rows, 0);
        // Every node receives only its estimation observations. Its value can
        // be used when truncating the same construction tree at a shallower K.
        for (int i = 0; i < ne; ++i) {
            int id = 0;
            while (true) {
                Node& node = builder.nodes[id];
                node.est.add(REAL(ey)[i], INTEGER(ed)[i], builder.xi);
                if (node.left < 0) break;
                id = REAL(ex)[i + ne*node.feature] <= node.cut ? node.left : node.right;
            }
        }
        const int m = static_cast<int>(builder.nodes.size());
        SEXP result = PROTECT(Rf_allocMatrix(REALSXP, m, 13));
        for (int i = 0; i < m; ++i) {
            Node& node = builder.nodes[i];
            node.value = estimate(node.est, builder.rule);
            if (!std::isfinite(node.value)) {
                UNPROTECT(1);
                throw std::runtime_error("Nonfinite leaf estimate");
            }
            double row[] = {static_cast<double>(node.left + 1), static_cast<double>(node.right + 1),
                static_cast<double>(node.feature + 1), node.cut, static_cast<double>(node.depth),
                node.value, static_cast<double>(node.train.n), static_cast<double>(node.train.n1),
                static_cast<double>(node.est.n), static_cast<double>(node.est.n1), node.gain,
                static_cast<double>(node.train.n - node.train.n1),
                static_cast<double>(node.est.n - node.est.n1)};
            for (int j = 0; j < 13; ++j) REAL(result)[i + m*j] = row[j];
        }
        UNPROTECT(1);
        return result;
    } catch (const std::exception& e) { Rf_error("cky_fit: %s", e.what()); }
    return R_NilValue;
}

extern "C" SEXP cky_predict(SEXP tree, SEXP x, SEXP depths_) {
    const int m = Rf_nrows(tree), n = Rf_nrows(x), nd = Rf_length(depths_);
    const double *t = REAL(tree), *xx = REAL(x);
    const int *depths = INTEGER(depths_);
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, n, nd));
    for (int i = 0; i < n; ++i) {
        int id = 0;
        for (int k = 0; k < nd; ++k) {
            while (t[id] > 0 && t[id + m*4] < depths[k]) {
                const int j = static_cast<int>(t[id + m*2]) - 1;
                id = static_cast<int>(xx[i + n*j] <= t[id + m*3] ? t[id] : t[id + m]) - 1;
            }
            REAL(result)[i + n*k] = t[id + m*5];
        }
    }
    UNPROTECT(1);
    return result;
}

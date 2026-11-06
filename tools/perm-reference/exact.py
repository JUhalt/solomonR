# Exact and high-precision helpers for the perm_solomon() reference set
# (issue #134). Python standard library only; run as
#   python -I exact.py <mode> <input.json> <output.json>
# by lib.R. Nothing here knows the package or its tie rule: values are
# compared exactly (rationals, integer relations) or at 50 significant
# digits, never with a tolerance that decides a tie.
#
# Numbers arrive as strings "m/d" (integers m and d; for data d is a power of
# two, so every stored double is represented exactly) and leave as strings,
# codes, or relations (-1, 0, 1) to an observed class.
#
# Modes
#   gauss    Gaussian identity model, each class in exact rationals: L'b, its
#            HC3 variance, codes (NE, H1, ZZ, INF-C, ZERO, finite), and the
#            dense rank of |T| over the classes, per contrast and statistic.
#   struct   Rank facts of a design: whether each contrast is estimable, and
#            leverage 1 (deleting a row lowers the rank and the contrast is
#            then not estimable).
#   mlecert  Certificate that the maximum-likelihood estimate exists: a mean
#            vector strictly inside (0, 1) or (0, Inf) with X'F mu = X'F y,
#            built exactly from the fitted means of a converged fit.
#   dircert  Certificate of a direction of recession: an exact rational d with
#            x_i'd = 0 on the rows kept and the strict sign on the rows that
#            diverge (Albert & Anderson, 1984; Geyer, 2009).
#   satcmp   Saturated cells of models fitted by iteration (binomial logit,
#            Poisson log, negative binomial log with theta fixed and no
#            exposure): each contrast is log R with R rational and its HC3
#            variance V rational. Equality is decided exactly; for the
#            studentized statistic by the integer relation R_a^q = R_b^(+-p)
#            when V_a / V_b = (p/q)^2, which is complete by the
#            Gelfond-Schneider theorem (Baker, 1975). Order from doubles when
#            far apart, otherwise at 60 digits.
#   band     Newton's method at 70 digits for classes whose order relative to
#            an observed class doubles could not certify; returns the
#            statistics to 50 significant digits.
#   cluster  Cluster-level references with exact Stage 1 scores (least
#            squares, or condition means and rates) and every allocation.

import sys
import json
import math
from fractions import Fraction as F
from decimal import Decimal as D, getcontext
from functools import cmp_to_key
from itertools import combinations, product

INF_RANK = 10 ** 9


class Unexplained(Exception):
    pass


def fr(s):
    if isinstance(s, int):
        return F(s)
    if isinstance(s, float):
        return F(s)
    if "/" in s:
        m, d = s.split("/")
        return F(int(m), int(d))
    return F(int(s))


def fstr(x):
    return "%d/%d" % (x.numerator, x.denominator)


# ---------------------------------------------------------------------------
# Exact linear algebra on lists of Fractions.

def rref(rows, ncol):
    """Reduced row echelon form of the row space: (basis rows, pivot columns)."""
    M = [list(r) for r in rows if any(v != 0 for v in r)]
    piv = []
    r = 0
    for c in range(ncol):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [v / pv for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return M[:r], piv


def in_rowspace(basis, piv, v):
    v = list(v)
    for row, c in zip(basis, piv):
        if v[c] != 0:
            f = v[c]
            v = [a - f * b for a, b in zip(v, row)]
    return all(x == 0 for x in v)


def nullspace(rows, ncol):
    basis, piv = rref(rows, ncol)
    free = [c for c in range(ncol) if c not in piv]
    out = []
    for fc in free:
        v = [F(0)] * ncol
        v[fc] = F(1)
        for row, pc in zip(basis, piv):
            v[pc] = -row[fc]
        out.append(v)
    return out


def inverse(A):
    n = len(A)
    M = [list(A[i]) + [F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


def matvec(A, v):
    return [sum(a * b for a, b in zip(row, v)) for row in A]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def dense_ranks(n, items, cmp):
    ranks = [None] * n
    finite = []
    for i, key in items:
        if key is None:
            continue
        if key == "inf":
            ranks[i] = INF_RANK
        elif key == 0:
            ranks[i] = 0
        else:
            finite.append((i, key))
    finite.sort(key=cmp_to_key(lambda a, b: cmp(a[1], b[1])))
    r = 0
    prev = None
    for i, key in finite:
        if prev is None or cmp(prev, key) != 0:
            r += 1
        ranks[i] = r
        prev = key
    return ranks


def frac_cmp(a, b):
    return (a > b) - (a < b)


# ---------------------------------------------------------------------------
# gauss

def gauss(inp):
    contrasts = {k: [fr(v) for v in L] for k, L in inp["contrasts"].items()}
    names = list(contrasts.keys())
    k = inp["k"]
    ncol = 4 + k
    classes = inp["classes"]
    n = len(classes)
    out = {c: {"diff": {"code": [], "approx": [], "key": []},
               "stud": {"code": [], "approx": [], "key": []}} for c in names}
    for cl in classes:
        X = []
        for cell, z in zip(cl["cell"], cl["z"]):
            row = [F(0)] * 4
            row[cell] = F(1)
            row += [fr(v) for v in z]
            X.append(row)
        y = [fr(v) for v in cl["y"]]
        f = [F(int(v)) for v in cl["f"]]
        basis, piv = rref(X, ncol)
        Xb = [[row[j] for j in piv] for row in X]
        r = len(piv)
        m = len(Xb)
        A = [[sum(f[i] * Xb[i][a] * Xb[i][b] for i in range(m)) for b in range(r)] for a in range(r)]
        Ai = inverse(A)
        Xty = [sum(f[i] * Xb[i][a] * y[i] for i in range(m)) for a in range(r)]
        beta = matvec(Ai, Xty)
        e = [yi - dot(xb, beta) for yi, xb in zip(y, Xb)]
        AiX = [matvec(Ai, xb) for xb in Xb]
        h = [dot(xb, ax) for xb, ax in zip(Xb, AiX)]
        for c in names:
            L = contrasts[c]
            if not in_rowspace(basis, piv, L):
                # NE: a cell the contrast uses has no rows; ALIAS: a
                # covariate is collinear with the cells under this labeling
                empty = any(L[j] != 0 and all(row[j] == 0 for row in X) for j in range(4))
                for s in ("diff", "stud"):
                    out[c][s]["code"].append("NE" if empty else "ALIAS")
                    out[c][s]["approx"].append(None)
                    out[c][s]["key"].append(None)
                continue
            Lb = [L[j] for j in piv]
            g = matvec(Ai, Lb)  # Ai is symmetric
            est = dot(g, Xty)
            V = F(0)
            h1 = False
            for xb, ei, hi, fi in zip(Xb, e, h, f):
                ci = dot(g, xb)
                if hi == 1:
                    if ci != 0:
                        h1 = True
                    continue
                t = ci * ei / (1 - hi)
                V += fi * t * t
            if est == 0:
                out[c]["diff"]["code"].append("ZERO")
                out[c]["diff"]["key"].append(0)
            else:
                out[c]["diff"]["code"].append("finite")
                out[c]["diff"]["key"].append(abs(est))
            out[c]["diff"]["approx"].append(float(est))
            if h1:
                code, key, ap = "H1", None, None
            elif V == 0:
                if est == 0:
                    code, key, ap = "ZZ", 0, 0.0
                else:
                    code, key, ap = "INF-C", "inf", math.copysign(math.inf, float(est))
            elif est == 0:
                code, key, ap = "ZERO", 0, 0.0
            else:
                z2 = est * est / V
                code, key = "finite", z2
                ap = math.copysign(math.sqrt(float(z2)), float(est))
            out[c]["stud"]["code"].append(code)
            out[c]["stud"]["key"].append(key)
            out[c]["stud"]["approx"].append(None if ap is None or math.isinf(ap) else ap)
    res = {}
    for c in names:
        res[c] = {}
        for s in ("diff", "stud"):
            keys = out[c][s]["key"]
            ranks = dense_ranks(n, list(enumerate(keys)), frac_cmp)
            res[c][s] = {"code": out[c][s]["code"], "approx": out[c][s]["approx"], "rank": ranks}
    return res


# ---------------------------------------------------------------------------
# struct

def struct(inp):
    out = []
    for u in inp["units"]:
        X = [[fr(v) for v in row] for row in u["X"]]
        ncol = len(X[0]) if X else 0
        basis, piv = rref(X, ncol)
        rX = len(piv)
        drop_rank = []
        for i in range(len(X)):
            if u.get("f") is not None and int(u["f"][i]) > 1:
                continue  # a row type with copies cannot have leverage 1
            Xi = X[:i] + X[i + 1:]
            bi, pi = rref(Xi, ncol)
            if len(pi) < rX:
                drop_rank.append((i, bi, pi))
        res = {"rank": rX, "piv": piv, "c": {}}
        for name, Ls in u["L"].items():
            L = [fr(v) for v in Ls]
            est = in_rowspace(basis, piv, L)
            h1 = False
            if est:
                for i, bi, pi in drop_rank:
                    if not in_rowspace(bi, pi, L):
                        h1 = True
                        break
            ncell = u.get("ncell", 4)
            alias = (not est) and not any(L[j] != 0 and all(row[j] == 0 for row in X)
                                          for j in range(min(ncell, len(L))))
            res["c"][name] = {"est": est, "h1": h1, "alias": alias}
        res["lev1_rows"] = [i for i, _, _ in drop_rank]
        out.append(res)
    return out


# ---------------------------------------------------------------------------
# mlecert and dircert

def mlecert_one(it):
    fam = it["fam"]
    X = [[fr(v) for v in row] for row in it["X"]]
    y = [fr(v) for v in it["y"]]
    f = [F(int(v)) for v in it["f"]]
    mu = [fr(v) for v in it["mu"]]
    n = len(X)
    ncol = len(X[0])
    basis, piv = rref(X, ncol)
    Xb = [[row[j] for j in piv] for row in X]
    r = len(piv)
    # mu2 = mu + D X c with D = diag(var(mu)) > 0 and X'F mu2 = X'F y: a
    # multiplicative correction, which keeps mu2 inside the bounds when the
    # fitted means are close to the estimate
    dvar = [m * (1 - m) if fam == "binomial" else m for m in mu]
    A = [[sum(f[i] * dvar[i] * Xb[i][a] * Xb[i][b] for i in range(n)) for b in range(r)] for a in range(r)]
    rhs = [sum(f[i] * Xb[i][a] * (y[i] - mu[i]) for i in range(n)) for a in range(r)]
    c = matvec(inverse(A), rhs)
    mu2 = [mu[i] + dvar[i] * dot(Xb[i], c) for i in range(n)]
    ok = all(m > 0 for m in mu2)
    if fam == "binomial":
        ok = ok and all(m < 1 for m in mu2)
    # exact check of the moment equations
    for a in range(r):
        if sum(f[i] * Xb[i][a] * (mu2[i] - y[i]) for i in range(n)) != 0:
            ok = False
    return ok


def mlecert(inp):
    return [mlecert_one(it) for it in inp["items"]]


def dircert_one(it):
    fam = it["fam"]
    X = [[fr(v) for v in row] for row in it["X"]]
    y = [fr(v) for v in it["y"]]
    I0 = set(int(i) for i in it["I0"])
    d0 = [float(v) for v in it["d"]]
    ncol = len(X[0])
    keep = [X[i] for i in range(len(X)) if i not in I0]
    N = nullspace(keep, ncol) if keep else [[F(1) if j == k else F(0) for j in range(ncol)]
                                              for k in range(ncol)]
    if not N:
        return False
    # least squares coefficients of d0 on the basis, in floats, then rounded
    Nf = [[float(v) for v in col] for col in N]
    k = len(Nf)
    G = [[sum(Nf[a][j] * Nf[b][j] for j in range(ncol)) for b in range(k)] for a in range(k)]
    rhs = [sum(Nf[a][j] * d0[j] for j in range(ncol)) for a in range(k)]
    # solve G c = rhs (small, by Gaussian elimination in floats)
    M = [G[i] + [rhs[i]] for i in range(k)]
    for c in range(k):
        p = max(range(c, k), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        if M[c][c] == 0:
            return False
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(k):
            if r != c:
                fct = M[r][c]
                M[r] = [a - fct * b for a, b in zip(M[r], M[c])]
    coef = [M[i][k] for i in range(k)]
    for den in (10 ** 3, 10 ** 6, 10 ** 9, 10 ** 12):
        cq = [F(v).limit_denominator(den) for v in coef]
        d = [sum(cq[a] * N[a][j] for a in range(k)) for j in range(ncol)]
        good = True
        for i in I0:
            v = dot(X[i], d)
            if fam == "binomial":
                if (y[i] == 1 and not v > 0) or (y[i] == 0 and not v < 0):
                    good = False
                    break
            else:
                if y[i] != 0 or not v < 0:
                    good = False
                    break
        if good:
            return True
    return False


def dircert(inp):
    return [dircert_one(it) for it in inp["items"]]


# ---------------------------------------------------------------------------
# satcmp: saturated cells.

def arm_parts(fam, rows):
    """rows: list of (y, e, f). Returns (rate or odds, HC3 variance of its log)."""
    ys = [fr(r[0]) for r in rows]
    es = [fr(r[1]) for r in rows]
    fs = [int(r[2]) for r in rows]
    n = sum(fs)
    S = sum(y * f for y, f in zip(ys, fs))
    if fam == "binomial":
        p = S / n
        odds = p / (1 - p)
        V = 1 / (n * p * (1 - p) * (1 - F(1, n)) ** 2)
        return odds, V
    E = sum(e * f for e, f in zip(es, fs))
    rate = S / E
    # Poisson log, and the negative binomial log with theta fixed and no
    # exposure (theta cancels): bread 1/S, leverage e_i/E.
    V = F(0)
    for y, e, f in zip(ys, es, fs):
        V += f * ((y - e * rate) / (1 - e / E)) ** 2
    V = V / (S * S)
    return rate, V


def dec(x):
    return D(x.numerator) / D(x.denominator)


def stud_equal(a, b):
    Ra, Va = a
    Rb, Vb = b
    if Ra == 1 or Rb == 1:
        return Ra == 1 and Rb == 1
    q2 = Va / Vb
    num, den = q2.numerator, q2.denominator
    rn, rd = math.isqrt(num), math.isqrt(den)
    if rn * rn != num or rd * rd != den:
        return False
    # |log Ra| / sqrt(Va) = |log Rb| / sqrt(Vb)  <=>  |log Ra| = (rn/rd) |log Rb|
    #   <=>  Ra^rd = Rb^(+-rn)
    if max(rn, rd) > 4000:
        getcontext().prec = 120
        la = abs(dec(Ra).ln()) * D(rd)
        lb = abs(dec(Rb).ln()) * D(rn)
        getcontext().prec = 60
        if abs(la - lb) <= D(10) ** -100 * max(la, lb):
            raise Unexplained("integer relation with large exponents")
        return False
    A = Ra ** rd
    B = Rb ** rn
    return A == B or A * B == 1


def stud_cmp(a, b):
    """Compare |log R| / sqrt(V); a, b are (R, V) with R != 1 and V > 0."""
    if stud_equal(a, b):
        return 0
    Ra, Va = a
    Rb, Vb = b
    da = math.log(float(Ra)) ** 2 / float(Va)
    db = math.log(float(Rb)) ** 2 / float(Vb)
    if math.isfinite(da) and math.isfinite(db) and abs(da - db) >= 1e-9 * max(da, db):
        return (da > db) - (da < db)
    getcontext().prec = 60
    ta = dec(Ra).ln() ** 2 / dec(Va)
    tb = dec(Rb).ln() ** 2 / dec(Vb)
    if abs(ta - tb) <= D(10) ** -35 * max(ta, tb):
        raise Unexplained("saturated classes equal to 35 digits without an integer relation")
    return (ta > tb) - (ta < tb)


def diff_key(R):
    return R if R >= 1 else 1 / R


def satcmp(inp):
    fam = inp["fam"]
    parts = []
    for p in inp["parts"]:
        if p is None:
            parts.append(None)
            continue
        rt, vt = arm_parts(fam, p["t"])
        rc, vc = arm_parts(fam, p["c"])
        parts.append((rt / rc, vt + vc))
    res = {}
    obs = inp["obs"]
    for cname, signs in inp["contrasts"].items():
        vals = []
        for cl in inp["classes"]:
            R = F(1)
            V = F(0)
            for k, s in enumerate(signs):
                if s == 0:
                    continue
                pr = parts[cl[k]]
                R = R * pr[0] if s > 0 else R / pr[0]
                V += pr[1]
            vals.append((R, V))
        dcode = ["ZERO" if R == 1 else "finite" for R, V in vals]
        scode = []
        for R, V in vals:
            if V == 0:
                scode.append("ZZ" if R == 1 else "INF-C")
            else:
                scode.append("ZERO" if R == 1 else "finite")
        cache = {}
        drel = []
        srel = []
        for o in obs:
            Ro, Vo = vals[o]
            dr = []
            sr = []
            for i, (R, V) in enumerate(vals):
                # difference: |log R| against |log Ro|
                ki = diff_key(R)
                ko = diff_key(Ro)
                dr.append((ki > ko) - (ki < ko))
                # studentized
                ci, co = scode[i], scode[o]
                rank = {"ZERO": 0, "ZZ": 0, "INF-C": 2, "finite": 1}
                if ci != "finite" or co != "finite":
                    a, b = rank[ci], rank[co]
                    sr.append((a > b) - (a < b))
                    continue
                key = (R, V, Ro, Vo)
                if key not in cache:
                    cache[key] = stud_cmp((R, V), (Ro, Vo))
                sr.append(cache[key])
            drel.append(dr)
            srel.append(sr)
        res[cname] = {"diff": {"code": dcode, "rel": drel,
                               "approx": [math.log(float(R)) for R, V in vals]},
                      "stud": {"code": scode, "rel": srel,
                               "approx": [math.log(float(R)) / math.sqrt(float(V)) if V > 0 else None
                                          for R, V in vals]},
                      "R": [fstr(R) for R, V in vals], "V": [fstr(V) for R, V in vals]}
    return res


# ---------------------------------------------------------------------------
# band: Newton's method at 70 digits.

def dsolve(A, b):
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * bb for a, bb in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]


def dinverse(A):
    n = len(A)
    cols = []
    for j in range(n):
        e = [D(1) if i == j else D(0) for i in range(n)]
        cols.append(dsolve(A, e))
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def band_one(it):
    getcontext().prec = 70
    fam = it["fam"]
    X = [[dec(fr(v)) for v in row] for row in it["X"]]
    y = [dec(fr(v)) for v in it["y"]]
    f = [D(int(v)) for v in it["f"]]
    n = len(y)
    p = len(X[0])
    off = [D(0)] * n
    if it.get("e") is not None:
        off = [dec(fr(v)).ln() for v in it["e"]]
    theta = dec(fr(it["theta"])) if it.get("theta") is not None else None
    b = [dec(fr(v)) for v in it["b0"]]

    def parts(b):
        eta = [sum(xij * bj for xij, bj in zip(xi, b)) + oi for xi, oi in zip(X, off)]
        if fam == "gaussian":
            mu = eta
            w = [D(1)] * n
            s = [yi - m for yi, m in zip(y, mu)]
            hw = w
        elif fam == "binomial":
            mu = [1 / (1 + (-e).exp()) for e in eta]
            w = [m * (1 - m) for m in mu]
            s = [yi - m for yi, m in zip(y, mu)]
            hw = w
        elif fam == "poisson":
            mu = [e.exp() for e in eta]
            w = mu
            s = [yi - m for yi, m in zip(y, mu)]
            hw = w
        else:
            mu = [e.exp() for e in eta]
            w = [m * theta / (theta + m) for m in mu]
            s = [(yi - m) * theta / (theta + m) for yi, m in zip(y, mu)]
            hw = [(yi + theta) * theta * m / (theta + m) ** 2 for yi, m in zip(y, mu)]
        return mu, w, s, hw

    size = None
    for _ in range(80):
        mu, w, s, hw = parts(b)
        score = [sum(f[i] * X[i][j] * s[i] for i in range(n)) for j in range(p)]
        H = [[sum(f[i] * X[i][a] * hw[i] * X[i][c] for i in range(n)) for c in range(p)]
             for a in range(p)]
        step = dsolve(H, score)
        b = [bj + sj for bj, sj in zip(b, step)]
        size = max(abs(v) for v in step) / max(D(1), max(abs(v) for v in b))
        if size < D(10) ** -60:
            break
    if size is None or size > D(10) ** -55:
        raise ValueError("70-digit Newton did not converge: %s" % size)
    mu, w, s, hw = parts(b)
    A = [[sum(f[i] * X[i][a] * w[i] * X[i][c] for i in range(n)) for c in range(p)] for a in range(p)]
    Ai = dinverse(A)
    AiX = [[sum(Ai[a][c] * X[i][c] for c in range(p)) for a in range(p)] for i in range(n)]
    h = [w[i] * sum(X[i][a] * AiX[i][a] for a in range(p)) for i in range(n)]
    out = {}
    skip = set(int(v) for v in (it.get("skip") or []))
    for name, spec in it["parts"].items():
        L = [dec(fr(v)) for v in spec["L"]]
        wown = dec(fr(spec["w"]))
        est = sum(lj * bj for lj, bj in zip(L, b))
        V = D(0)
        for i in range(n):
            if i in skip:
                continue  # leverage exactly 1 and c_i exactly 0: adds 0
            ci = sum(L[a] * AiX[i][a] for a in range(p))
            t = ci * s[i] / (1 - h[i])
            V += f[i] * t * t
        num = wown * est
        var = wown * wown * V
        for fx in spec.get("fixed", []):
            wk = dec(fr(fx["w"]))
            num += wk * dec(fr(fx["R"])).ln()
            var += wk * wk * dec(fr(fx["V"]))
        T = num / var.sqrt() if var > 0 else None
        out[name] = {"diff": format(+num, ".50e"),
                     "stud": format(+T, ".50e") if T is not None else None}
    return out


def band(inp):
    out = []
    for it in inp["items"]:
        try:
            out.append(band_one(it))
        except Exception as e:  # reported; the classes stay undecided
            out.append({"error": repr(e)})
    return out


# ---------------------------------------------------------------------------
# cluster: exact Stage 1 scores and allocation statistics.

def cluster(inp):
    rows = inp["rows"]          # list of {cl, y, p, z: [...], e}
    stage1 = inp["stage1"]      # "ols" or "means"
    design = inp["design"]      # "A" or "B"
    ys = [fr(r["y"]) for r in rows]
    es = [fr(r["e"]) for r in rows]
    pf = [r["p"] for r in rows]
    n = len(rows)
    if stage1 == "ols":
        X = []
        for r in rows:
            row = [F(1) if r["p"] == 0 else F(0), F(1) if r["p"] == 1 else F(0)]
            row += [fr(v) for v in r["z"]]
            X.append(row)
        ncol = len(X[0])
        basis, piv = rref(X, ncol)
        Xb = [[row[j] for j in piv] for row in X]
        rk = len(piv)
        A = [[sum(Xb[i][a] * Xb[i][b] for i in range(n)) for b in range(rk)] for a in range(rk)]
        Ai = inverse(A)
        beta = matvec(Ai, [sum(Xb[i][a] * ys[i] for i in range(n)) for a in range(rk)])
        fitted = [dot(xb, beta) for xb in Xb]
    else:
        tot = {}
        for y, e, p in zip(ys, es, pf):
            S, E = tot.get(p, (F(0), F(0)))
            tot[p] = (S + y, E + e)
        fitted = []
        for e, p in zip(es, pf):
            S, E = tot[p]
            fitted.append(e * S / E)
    resid = [y - f for y, f in zip(ys, fitted)]
    sums = {}
    for r, ei, e in zip(rows, resid, es):
        key = (r["cl"], r["p"]) if design == "B" else r["cl"]
        s, z = sums.get(key, (F(0), F(0)))
        sums[key] = (s + ei, z + e)
    score = {k: v[0] / v[1] for k, v in sums.items()}
    strata = []
    for st in inp["strata"]:
        ids = st["ids"]
        if design == "B":
            wP, wU = fr(inp["wB"][0]), fr(inp["wB"][1])
            s = [score[(c, 1)] * wP + score[(c, 0)] * wU for c in ids]
        else:
            s = [score[c] for c in ids]
        strata.append({"ids": ids, "s": s, "t": st["t"], "w": fr(st["w"])})
    per = [list(combinations(range(len(st["ids"])), sum(st["t"]))) for st in strata]
    keys = {"diff": [], "stud": []}
    codes = {"diff": [], "stud": []}
    approx = {"diff": [], "stud": []}
    alloc_out = []
    obs_index = None
    for a in product(*[range(len(p)) for p in per]):
        est = F(0)
        var = F(0)
        tsets = []
        for st, pc, j in zip(strata, per, a):
            T = set(pc[j])
            tsets.append([st["ids"][i] for i in sorted(T)])
            st_t = [st["s"][i] for i in range(len(st["s"])) if i in T]
            st_c = [st["s"][i] for i in range(len(st["s"])) if i not in T]
            nt, nc = len(st_t), len(st_c)
            mt = sum(st_t) / nt
            mc = sum(st_c) / nc
            vt = sum((v - mt) ** 2 for v in st_t) / (nt - 1) if nt > 1 else None
            vc = sum((v - mc) ** 2 for v in st_c) / (nc - 1) if nc > 1 else None
            est += st["w"] * (mt - mc)
            if vt is None or vc is None:
                var = None
            elif var is not None:
                var += st["w"] ** 2 * (vt / nt + vc / nc)
        alloc_out.append(tsets)
        if all(sorted(st["ids"][i] for i in range(len(st["t"])) if st["t"][i] == 1) == ts
               for st, ts in zip(strata, tsets)):
            obs_index = len(alloc_out) - 1
        if est == 0:
            codes["diff"].append("ZERO"); keys["diff"].append(0)
        else:
            codes["diff"].append("finite"); keys["diff"].append(abs(est))
        approx["diff"].append(float(est))
        if var is None:
            codes["stud"].append("NA"); keys["stud"].append(None); approx["stud"].append(None)
        elif var == 0:
            if est == 0:
                codes["stud"].append("ZZ"); keys["stud"].append(0); approx["stud"].append(0.0)
            else:
                codes["stud"].append("INF-C"); keys["stud"].append("inf"); approx["stud"].append(None)
        elif est == 0:
            codes["stud"].append("ZERO"); keys["stud"].append(0); approx["stud"].append(0.0)
        else:
            z2 = est * est / var
            codes["stud"].append("finite"); keys["stud"].append(z2)
            approx["stud"].append(math.copysign(math.sqrt(float(z2)), float(est)))
    res = {"alloc": alloc_out, "obs": obs_index,
           "scores": {str(k): fstr(v) for k, v in score.items()}}
    for s in ("diff", "stud"):
        res[s] = {"code": codes[s], "approx": approx[s],
                  "rank": dense_ranks(len(alloc_out), list(enumerate(keys[s])), frac_cmp)}
    return res


# ---------------------------------------------------------------------------
# satparts: (R, V) of single saturated conditions, as exact rationals.

def satparts(inp):
    fam = inp["fam"]
    out = []
    for p in inp["parts"]:
        rt, vt = arm_parts(fam, p["t"])
        rc, vc = arm_parts(fam, p["c"])
        out.append({"R": fstr(rt / rc), "V": fstr(vt + vc)})
    return out


# decgap: sign of |a| - |b| and their relative gap, for 50-digit strings.

def decgap(inp):
    getcontext().prec = 80
    out = []
    for a, b in inp["pairs"]:
        x = abs(D(a))
        y = abs(D(b))
        if x == y:
            out.append({"sign": 0, "rel": 0.0})
            continue
        rel = abs(x - y) / max(x, y)
        out.append({"sign": 1 if x > y else -1, "rel": float(rel)})
    return out


# bigsum: exact sums of products of integers given as decimal strings.

def bigsum(inp):
    out = []
    for it in inp["items"]:
        tot = sum(int(a) * int(b) for a, b in zip(it["a"], it["b"]))
        out.append(str(tot))
    return out


MODES = {"gauss": gauss, "struct": struct, "mlecert": mlecert, "dircert": dircert,
         "satcmp": satcmp, "satparts": satparts, "band": band, "cluster": cluster,
         "bigsum": bigsum, "decgap": decgap}


def run_one(mode, fin, fout):
    with open(fin) as f:
        inp = json.load(f)
    try:
        res = {"ok": True, "result": MODES[mode](inp)}
    except Unexplained as e:
        res = {"ok": False, "unexplained": str(e)}
    with open(fout, "w") as f:
        json.dump(res, f)


def main():
    if sys.argv[1] == "serve":
        # one request per line on stdin: "<mode> <input.json> <output.json>";
        # the reply "done" (or "fail <message>") is written once the output
        # file is complete
        for line in sys.stdin:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                break
            try:
                run_one(*parts)
                sys.stdout.write("done\n")
            except Exception as e:  # reported to the caller, which stops
                sys.stdout.write("fail %s\n" % repr(e).replace("\n", " "))
            sys.stdout.flush()
        return
    run_one(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()

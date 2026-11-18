// ============================================================================
//  nospec_twosteps.cpp — no q_spec (two-step method): τ search + reuse of SCRM's public two-step API
//  See nospec_twosteps.h and docs/NoSpec_no_specular_peak_two_step_design_algorithm_wiki.md
//
//  This file is a **fully standalone new implementation**: it modifies none of scrm.cpp / nospec.cpp / specindex*
//  or any other existing algorithm file, calling them only through their public interfaces. The two hard
//  criteria in SCRM for "q_s is observed" are disabled by setting scrm::Params::dqspec to a huge value (that
//  field is already a public parameter), with no switch added to SCRM.
//  Parameters: nospec2::Params is this mode's own complete parameter set (visible and tunable item by item in the UI); this file only
//  performs a single call adaptation in makeScrmParams, embedding no scrm::Params and reusing none of its default search ranges.
// ============================================================================
#include "nospec_twosteps.h"
#include "scrm.h" // Used only in this file: makeScrmParams below performs a single call adaptation

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_set>
#include <vector>

namespace gid {
namespace nospec2 {
namespace {

constexpr double kBig = 1e30;
constexpr double kPhi = 0.6180339887498949;

struct Range {
    double lo = 0.0, hi = 0.0;
};

// Search interval for τ derived from the peak table. q_s is the length of the contact-plane reciprocal vector, usually of the same order as q_z,
// but it may exceed all measured q_z (e.g. DIP: q_s=1.773 > max(q_z)=1.232), so the upper bound is
// 2.5·max(q_z); then 1.2·max|q| covers the case where "all peak q_z are very small".
Range autoRange(const std::vector<Peak>& pk) {
    double qzMax = 0.0, qMax = 0.0, qzMin = kBig;
    for (const Peak& p : pk) {
        qMax = std::max(qMax, std::hypot(p.qxy, p.qz));
        qzMax = std::max(qzMax, p.qz);
        if (p.qz > 1e-9) qzMin = std::min(qzMin, p.qz);
    }
    if (qzMin >= kBig * 0.5) qzMin = qMax;
    double lo = std::max(0.5 * qzMin, 0.15 * qMax);
    double hi = std::max(2.5 * qzMax, 1.2 * qMax);
    if (!(lo > 1e-6)) lo = 0.1;
    if (!(hi > lo * 1.01)) hi = lo * 4.0;
    return Range{lo, hi};
}

// τ scoring: for a given τ, (22) is a 3×3 linear system; take the minimum rmsQxy of the SCRM in-plane candidates.
// This step is cheap but has limited discriminative power—step one can find some (u,v,h,k) assignment for any τ
// that fits q_xy, and an overfitting valley at a wrong τ (DIP measured τ≈0.90 gives 0.0038) may beat the true value
// (τ≈1.773 gives 0.0054). Therefore the outer loop does two more things:
//   (1) keep multiple local minima and fill in τ candidates by scale diversity, not just the few smallest R;
//   (2) if the peak table contains rod peaks with q_xy≈0, send q_z,rod/m to step two as soft prior candidates.
// The final adjudication is always done by step two's integer consistency and global rmsQxyz (design draft §4.5.1/§4.5.2).
double tauScore(const std::vector<Peak>& pk, const scrm::Params& base, double tau,
                const Progress& quiet, Stats* st) {
    if (!(tau > 0.0)) return kBig;
    scrm::Params p = base;
    p.qspec = tau;
    scrm::Stats ips;
    const std::vector<scrm::InPlaneCand> cands = scrm::runInPlane(pk, p, quiet, &ips);
    if (st) st->ipCands += (long long)cands.size();
    double best = kBig;
    for (const scrm::InPlaneCand& c : cands) best = std::min(best, c.rmsQxy);
    return best;
}

// Solution deduplication key: cell parameters quantised to 0.02 Å / 0.1°, plus contact plane and q_spec (0.005 Å⁻¹).
// Unlike SCRM's cellKeyOf: q_spec must be included here, otherwise the same unit cell under different τ
// branches would be treated as different solutions and output repeatedly.
std::string solutionKey(const Solution& s) {
    char buf[256];
    std::snprintf(buf, sizeof(buf), "%lld/%lld/%lld/%lld/%lld/%lld/%d/%d/%d/%lld",
                  (long long)std::llround(s.cell.a * 50.0),
                  (long long)std::llround(s.cell.b * 50.0),
                  (long long)std::llround(s.cell.c * 50.0),
                  (long long)std::llround(s.cell.al * 10.0),
                  (long long)std::llround(s.cell.be * 10.0),
                  (long long)std::llround(s.cell.ga * 10.0), s.uvw.h, s.uvw.k, s.uvw.l,
                  (long long)std::llround(s.qspec * 200.0));
    return std::string(buf);
}

// The only "parameter adaptation" entry: maps this mode's own Params to the call parameters of the SCRM two-step API.
// Deliberately does not rely on scrm::Params default-constructed values to fill any tunable item—all search ranges, criteria, and
// reduction switches are explicitly given by nospec2::Params (visible item by item in the UI), and the mapping only moves fields. Three fixed semantics:
//   qspec = current τ; dqspec = kBig (no measured q_s, disables the hard self-consistency at scrm.cpp:834/:1143);
//   earlyStop = false (multiple τ run together; early exit on any τ would break the global ordering).
// This is also the only contact point between this file and scrm::Params; neither SCRM source nor its default configuration is modified.
scrm::Params makeScrmParams(const Params& q, double tau) {
    scrm::Params p;
    p.qspec = tau;
    p.system = q.system;
    p.fixedUvw = false;  // contact plane (u,v,w) is given by enumeration, not fixed
    p.fu = 0; p.fv = 0; p.fw = 1; // unused when fixedUvw=false; given explicitly, not relying on default construction
    p.uvLimit = q.uvLimit;
    p.hkLimit1 = q.hk1;
    p.hkLimit2 = q.hk2;
    p.seed3Half = q.seed3Half;
    p.linesToPermute = q.seedPeaks;
    p.redundant = q.redundant;
    p.maxPart1 = q.maxPart1;
    p.seedL = q.seedL;
    p.lLimit = q.lLimit;
    p.wLimit = q.wLimit;
    p.maxTrials = q.maxTrials;
    p.refineTop = q.refineTop;
    p.maxSolutions = q.maxSolutions;
    p.earlyStop = false;
    p.dqxyCutoff = q.dqxyCutoff;
    p.dqzCutoff = q.dqzCutoff;
    p.dqxyzCutoff = q.dqxyzCutoff;
    p.dqspec = kBig;
    p.restrictCell = q.restrictCell;
    p.aLo = q.aLo;   p.aHi = q.aHi;
    p.bLo = q.bLo;   p.bHi = q.bHi;
    p.cLo = q.cLo;   p.cHi = q.cHi;
    p.alLo = q.alLo; p.alHi = q.alHi;
    p.beLo = q.beLo; p.beHi = q.beHi;
    p.gaLo = q.gaLo; p.gaHi = q.gaHi;
    p.volLo = q.volLo; p.volHi = q.volHi;
    p.niggli = q.niggli;
    p.threads = q.threads;
    return p;
}

} // namespace

Params Params::sanitized() const {
    Params p = *this;
    // τ search
    p.qsSteps = std::clamp(p.qsSteps, 4, 64);
    p.qsRefine = std::clamp(p.qsRefine, 0, 24);
    p.maxTau = std::clamp(p.maxTau, 1, 24);
    if (p.rodQxyTol > 0.0) p.rodQxyTol = std::clamp(p.rodQxyTol, 1e-6, 0.5);
    p.rodOrders = std::clamp(p.rodOrders, 1, 8);
    p.rodPrior = std::clamp(p.rodPrior, 0.05, 1.0);
    if (p.qsLo > 0.0 && p.qsHi > 0.0 && !(p.qsHi > p.qsLo)) p.qsHi = p.qsLo * 2.0;
    p.coarseHk1 = std::clamp(p.coarseHk1, 1, 8);
    p.coarseHk2 = std::clamp(p.coarseHk2, 1, 12);
    p.coarseSeedPeaks = std::clamp(p.coarseSeedPeaks, 1, 12);
    p.coarseRedundant = std::clamp(p.coarseRedundant, 1, 8);
    p.coarseMaxPart1 = std::clamp(p.coarseMaxPart1, 1, 2000);
    // Step one: in-plane enumeration
    p.uvLimit = std::clamp(p.uvLimit, 0, 5);
    p.hk1 = std::clamp(p.hk1, 1, 8);
    p.hk2 = std::clamp(p.hk2, 1, 12);
    p.seed3Half = std::clamp(p.seed3Half, 0, 12);
    p.seedPeaks = std::clamp(p.seedPeaks, 1, 12);
    p.redundant = std::clamp(p.redundant, 1, 8);
    p.maxPart1 = std::clamp(p.maxPart1, 1, 2000);
    // Step two: out-of-plane enumeration
    p.seedL = std::clamp(p.seedL, 0, 6);
    p.seedLFallback = std::clamp(p.seedLFallback, 0, 6);
    p.lLimit = std::clamp(p.lLimit, 1, 20);
    p.wLimit = std::clamp(p.wLimit, 1, 12);
    p.maxTrials = std::max(0, p.maxTrials);
    p.refineTop = std::clamp(p.refineTop, 1, 2000);
    // Criteria / output
    if (!(p.dqxyCutoff > 0.0)) p.dqxyCutoff = 0.05;
    if (!(p.dqzCutoff > 0.0)) p.dqzCutoff = 0.05;
    if (!(p.dqxyzCutoff > 0.0)) p.dqxyzCutoff = 0.05;
    if (p.aHi <= p.aLo) p.aHi = p.aLo + 1.0;
    if (p.bHi <= p.bLo) p.bHi = p.bLo + 1.0;
    if (p.cHi <= p.cLo) p.cHi = p.cLo + 1.0;
    if (p.volHi <= p.volLo) p.volHi = p.volLo + 1.0;
    p.maxIpCands = std::clamp(p.maxIpCands, 1, 64);
    p.maxSolutions = std::clamp(p.maxSolutions, 1, 500);
    return p;
}

std::vector<Solution> runNoSpecTwoSteps(const std::vector<Peak>& peaksIn, const Params& prmIn,
                                        const Progress& prog, Stats* stOut) {
    const auto t0 = std::chrono::steady_clock::now();
    Stats st;
    st.peaksIn = (int)peaksIn.size();

    std::vector<Peak> pk = peaksIn;
    normalizeQxySigns(pk);
    // Use the same filtering rule as SCRM's two-step functions, so InPlaneCand::hk indices stay consistent with the peak table.
    std::vector<Peak> data;
    data.reserve(pk.size());
    for (const Peak& q : pk)
        if (q.qxy > 1e-6 && std::isfinite(q.qxy) && std::isfinite(q.qz)) data.push_back(q);
    st.peaksUsed = (int)data.size();
    // Without specular peaks, step one has 4 continuous unknowns; the practical lower bound is 4 peaks (design draft §3.4).
    if (data.size() < 4) return {};

    const Params prm = prmIn.sanitized();
    Progress quiet;
    quiet.cancelled = prog.cancelled; // only forward the cancellation signal to the inner layer; progress is reported by the outer layer

    // ---- τ coarse scan parameters: all come from the visible/tunable coarse* items in Params, no hidden constants ----
    // The coarse scan evaluates dozens of τ with a lightweight budget; the final stage reruns fully with the main parameters below.
    // The placeholder qspec=1 is overwritten by the current τ inside tauScore.
    scrm::Params base = makeScrmParams(prm, 1.0);
    base.threads = 1; // avoid repeatedly creating threads for each τ
    base.hkLimit1 = prm.coarseHk1;
    base.hkLimit2 = prm.coarseHk2;
    base.linesToPermute = prm.coarseSeedPeaks;
    base.redundant = prm.coarseRedundant;
    base.maxPart1 = prm.coarseMaxPart1;

    const Range r = autoRange(data);
    const double lo = (prm.qsLo > 0.0) ? prm.qsLo : r.lo;
    const double hi = (prm.qsHi > 0.0) ? prm.qsHi : r.hi;
    st.qsLoActual = lo;
    st.qsHiActual = hi;
    if (!(hi > lo) || prog.stop()) return {};

    // ---- Coarse scan ----
    const int n = prm.qsSteps;
    std::vector<double> grid((size_t)n), score((size_t)n, kBig);
    for (int j = 0; j < n; ++j) {
        if (prog.stop()) return {};
        grid[(size_t)j] = lo + (hi - lo) * (double)j / (double)(n - 1);
        score[(size_t)j] = tauScore(data, base, grid[(size_t)j], quiet, &st);
        ++st.tauTried;
        prog(0.05 + 0.60 * (double)j / (double)(n - 1), "NoSpec2 Step one: scanning q_spec");
    }

    // ---- Local minimum refinement (golden section) ----
    // Design draft §4.5.2: for N=4, det[M₄(τ)|−b]=0 may have multiple roots, all fitting 4 peaks exactly;
    // here each local minimum of the in-plane residual generates one τ candidate, handed to step two's integer-consistency adjudication.
    std::vector<std::pair<double, double>> scored; // (R, τ)
    int bestJ = 0;
    for (int j = 1; j < n; ++j)
        if (score[(size_t)j] < score[(size_t)bestJ]) bestJ = j;
    scored.emplace_back(score[(size_t)bestJ], grid[(size_t)bestJ]);

    // Rod peak soft prior: peaks with q_xy≈0 are excluded from fitting (same as SCRM), but q_z,rod/m
    // are natural candidate positions for q_s. Add them to the τ candidates with a score discount so they enter step two first;
    // these are only candidates, not labels saying "this is the measured q_spec", and no hard constraint is imposed (§8 A5 / §7.5.4).
    if (prm.rodQxyTol > 0.0) {
        for (const Peak& q : pk) {
            if (!(q.qz > 0.0) || !(std::fabs(q.qxy) < prm.rodQxyTol)) continue;
            for (int m = 1; m <= prm.rodOrders; ++m) {
                const double t = q.qz / (double)m;
                if (!(t > lo && t < hi)) continue;
                const double f = tauScore(data, base, t, quiet, &st);
                ++st.tauTried;
                // The first-order rod peak is the most likely q_s, with the smallest discount (highest priority); higher orders m weaken in turn.
                // In measurements, putting q_z/2 first pushed the true τ to the second candidate and lost the manual solution.
                const double prior = std::min(1.0, prm.rodPrior + 0.2 * (double)(m - 1));
                scored.emplace_back(f * prior, t);
            }
        }
    }

    auto eval = [&](double t) -> double {
        const double f = tauScore(data, base, t, quiet, &st);
        ++st.tauTried;
        return f;
    };
    for (int j = 1; j + 1 < n; ++j) {
        if (score[(size_t)j] >= kBig * 0.5) continue;
        if (!(score[(size_t)j] <= score[(size_t)j - 1] &&
              score[(size_t)j] <= score[(size_t)j + 1]))
            continue;
        double a = grid[(size_t)j - 1], b = grid[(size_t)j + 1];
        double x1 = b - kPhi * (b - a), x2 = a + kPhi * (b - a);
        double f1 = eval(x1), f2 = eval(x2);
        for (int it = 0; it < prm.qsRefine && (b - a) > 1e-6; ++it) {
            if (f1 < f2) {
                b = x2; x2 = x1; f2 = f1;
                x1 = b - kPhi * (b - a);
                f1 = eval(x1);
            } else {
                a = x1; x1 = x2; f1 = f2;
                x2 = a + kPhi * (b - a);
                f2 = eval(x2);
            }
        }
        const double tm = 0.5 * (a + b);
        const double fm = eval(tm);
        double tb = tm, fb = fm;
        if (f1 < fb) { tb = x1; fb = f1; }
        if (f2 < fb) { tb = x2; fb = f2; }
        scored.emplace_back(fb, tb);
    }

    std::sort(scored.begin(), scored.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    // Deduplicate first (after refinement several minima may converge to the same τ)
    std::vector<std::pair<double, double>> uniq;
    for (const auto& s : scored) {
        if (!(s.first < kBig * 0.5) || !(s.second > 0.0)) continue;
        bool dup = false;
        for (const auto& u : uniq)
            if (std::fabs(u.second - s.second) <= 1e-3 * std::max(1.0, s.second)) {
                dup = true; break;
            }
        if (!dup) uniq.push_back(s);
    }
    if (uniq.empty()) return {};

    // Taking only the top few by in-plane residual misses the true τ: overfitting (u,v) at a wrong τ often gives a smaller
    // rmsQxy (measured DIP: 0.0038 at τ≈0.90 beats 0.0054 at the true τ≈1.773).
    // Therefore half the slots are chosen by best R, and the other half greedily filled by distance from the selected τ, ensuring coverage of different scales.
    std::vector<double> taus;
    const size_t byScore = std::max<size_t>(1, (size_t)prm.maxTau / 2);
    for (const auto& u : uniq) {
        if (taus.size() >= byScore) break;
        taus.push_back(u.second);
    }
    while (taus.size() < (size_t)prm.maxTau) {
        size_t pick = uniq.size();
        double bestSep = -1.0;
        for (size_t i = 0; i < uniq.size(); ++i) {
            const double t = uniq[i].second;
            bool used = false;
            double sep = kBig;
            for (double x : taus) {
                const double d = std::fabs(t - x);
                if (d <= 1e-3 * std::max(1.0, t)) { used = true; break; }
                sep = std::min(sep, d / std::max(1.0, t));
            }
            if (used) continue;
            if (sep > bestSep) { bestSep = sep; pick = i; }
        }
        if (pick == uniq.size()) break;
        taus.push_back(uniq[pick].second);
    }
    if (taus.empty()) return {};
    st.tauKept = (int)taus.size();
    st.tauUsed = taus;
    st.bestTau = uniq.front().second;
    st.bestTauScore = uniq.front().first;

    // ---- Step two: full search only for the surviving τ ----
    // dqspec = kBig is the key to "no measured q_spec": SCRM at scrm.cpp:834 and :1143
    // uses |q_s(solution) − q_s(candidate)| ≤ dqspec for hard self-consistency, which only holds when q_s is an observed value.
    // Once disabled, w / l integer consistency (the 0.35 criterion inside runNormal) and the final rms govern.
    std::vector<Solution> all;
    std::vector<scrm::InPlaneCand> firstCands;
    scrm::Params firstP;
    std::vector<double> tauBest(taus.size(), kBig);

    auto secondPass = [&](const scrm::Params& p,
                          const std::vector<scrm::InPlaneCand>& cands, size_t ti) {
        for (const scrm::InPlaneCand& c : cands) {
            if (prog.stop()) break;
            scrm::Stats ns;
            std::vector<Solution> sols = scrm::runNormal(data, c, p, &ns);
            st.secondTrials += ns.secondTrials;
            st.rawSol += ns.rawSol;
            for (Solution& s : sols) {
                tauBest[ti] = std::min(tauBest[ti], s.rmsQxyz);
                if (!(s.rmsQxyz <= prm.dqxyzCutoff)) continue;
                s.dqspec = 0.0; // no measured q_spec, so the difference is undefined
                all.push_back(std::move(s));
            }
        }
    };

    for (size_t ti = 0; ti < taus.size(); ++ti) {
        if (prog.stop()) break;
        const double tau = taus[ti];
        scrm::Params p = makeScrmParams(prm, tau);
        scrm::Stats ips;
        std::vector<scrm::InPlaneCand> cands = scrm::runInPlane(data, p, quiet, &ips);
        st.ipCands += ips.ipCands;
        std::sort(cands.begin(), cands.end(),
                  [](const scrm::InPlaneCand& a, const scrm::InPlaneCand& b) {
                      return a.rmsQxy < b.rmsQxy;
                  });
        if ((int)cands.size() > prm.maxIpCands) cands.resize((size_t)prm.maxIpCands);
        if (ti == 0) { firstCands = cands; firstP = p; }
        prog(0.70 + 0.28 * (double)(ti + 1) / (double)taus.size(),
             "NoSpec2 Step two: solving out-of-plane components");
        secondPass(p, cands, ti);
    }

    // Fallback: when seed l radius 1 misses (the true solution's |l| can reach 2~3), widen it to 3 only on the best τ.
    // Same idea as SCRM-driven seedL fallback rounds, but confined to this module and affecting no existing path.
    const double good = 0.1 * prm.dqxyzCutoff;
    bool hasGood = false;
    for (const Solution& s : all)
        if (s.rmsQxyz <= good) { hasGood = true; break; }
    if (!hasGood && !firstCands.empty() && prm.seedLFallback > 0 &&
        firstP.seedL < prm.seedLFallback && !prog.stop()) {
        scrm::Params p = firstP;
        p.seedL = prm.seedLFallback;
        secondPass(p, firstCands, 0);
    }

    st.tauBestRms = tauBest;

    // ---- Sorting / deduplication ----
    // Within the same residual tier, smaller cells rank first (suppressing overfit supercells), but the tier width must adapt to the best solution:
    // a fixed 0.25·tolerance = 0.0125 would put the true solution at rms 0.0039 and a spurious solution at rms 0.0122
    // in the same tier, then sort by volume and let the spurious one win. Here take min(0.25·tolerance, 2·best rms).
    double bestRms = kBig;
    for (const Solution& s : all) bestRms = std::min(bestRms, s.rmsQxyz);
    const double floorRms = (bestRms < kBig * 0.5)
                                ? std::min(0.25 * prm.dqxyzCutoff,
                                           std::max(1e-4, 2.0 * bestRms))
                                : 0.25 * prm.dqxyzCutoff;
    std::sort(all.begin(), all.end(), [floorRms](const Solution& a, const Solution& b) {
        const double sa = std::max(a.rmsQxyz, floorRms), sb = std::max(b.rmsQxyz, floorRms);
        if (std::fabs(sa - sb) > 1e-9) return sa < sb;
        return a.vol < b.vol;
    });
    std::vector<Solution> out;
    std::unordered_set<std::string> seen;
    for (Solution& s : all) {
        canonicalizeUvw(s);
        if (!seen.insert(solutionKey(s)).second) continue;
        out.push_back(std::move(s));
        if ((int)out.size() >= prm.maxSolutions) break;
    }

    st.finalSol = (int)out.size();
    st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (stOut) *stOut = st;
    prog(1.0, "NoSpec2 complete");
    return out;
}

} // namespace nospec2
} // namespace gid

// ============================================================================
//  specindex.cpp — see specindex.h
//
//  Dispatcher: Improved / PeriodVoteDdl / Ransac / Ills (and Original) all forward to
//  improved::runSpecIndex (specindex_improved.cpp); Scrm forwards to the standalone
//  scrm::runSpecIndex (scrm.cpp). The original GIDInd part1/part2 port and the old SCRM path were removed.
// ============================================================================
#include "specindex.h"

#include "scrm.h"

namespace gid {

std::vector<Solution> runSpecIndex(const std::vector<Peak>& peaks, const SpecParams& prm,
                                   const Progress& prog) {
    if (!(prm.qspec > 0.0)) return {};
    if (prm.part1Alg == Part1Algorithm::Scrm)
        return scrm::runSpecIndex(peaks, scrm::fromSpecParams(prm), prog);
    return improved::runSpecIndex(peaks, prm, prog);
}

} // namespace gid

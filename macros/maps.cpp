/**
 * @file maps.cpp
 * @brief Utility functions for transforming and mapping tau-related quantities.
 * @details Provides a logit transform for raw discriminator scores and helper
 *   functions that map reconstructed and generator-level tau objects to
 *   decay-mode or genmatch integer codes suitable for histogramming in RDataFrame.
 * @author Izaak Neutelings (March 2026)
 */
#include <cmath>
#include "ROOT/RVec.hxx"

/// @brief Apply a logit transform to a vector of raw discriminator scores.
/// @details Computes logit(x) = 2*atanh(2*x-1), which maps [0,1] to [-inf,+inf]
///   and spreads out values near 0 and 1 for better histogram resolution.
///   The output is clipped to [-yclip,yclip], including under-/overflow.
///   Entries with the default undefined value (x < -0.5) are set to -99.
/// @param x     Input vector of raw scores in [0,1]
/// @param yclip Symmetric clip value so the output is [-yclip,+yclip] (default: 10).
/// @return Vector of logit-transformed scores, clipped to [-yclip, yclip],
///         or -99 for entries with undefined input (x < -0.5).
ROOT::RVecD logit(const ROOT::RVecD& x, double yclip=10) {
  ROOT::RVecD y(x.size());
  double xclip = 2/(1+std::exp(-yclip+0.005))-1; // 0.005 to avoid under/over-flow; 0.9999087492 for yclip=10
  for (size_t i=0; i<x.size(); ++i) {
    if (x[i]<-0.5) // default if undefined
      y[i] = -99.;
    else // clip x so y = logit(x) is between [-yclip,yclip]
      y[i] = 2*std::atanh(std::min(std::max(2*x[i]-1,-xclip),xclip));
  }
  return y;
}

/// @brief Map each reconstructed tau to an integer decay-mode/genmatch code for histogramming.
/// @details For genuine tau candidates (genmatch==5) the GenVisTau decay mode is forwarded.
///   For fake tau candidates, the genmatch category is mapped to a dedicated integer code.
///   Output codes:
///   | Code | Meaning |
///   |------|---------|
///   | -2   | Unrecognized genmatch (should never occur) |
///   | -1   | Real tau: no Tau–GenVisTau match |
///   |  0   | Real tau: OneProng0PiZero |
///   |  1   | Real tau: OneProng1PiZero |
///   |  2   | Real tau: OneProng2PiZero |
///   |  10  | Real tau: ThreeProng0PiZero |
///   |  11  | Real tau: ThreeProng1PiZero |
///   |  15  | Real tau: Other |
///   |  16  | Fake tau: prompt or tauonic electron (genmatch 1 or 3) |
///   |  17  | Fake tau: prompt or tauonic muon (genmatch 2 or 4) |
///   |  18  | Fake tau: q/g jet (genmatch 0) |
/// @param gm     `Tau_genPartFlav` — generator-level flavour/genmatch flag per reco tau.
/// @param genidx `Tau_genPartIdx`  — index into the GenVisTau collection for each reco tau.
/// @param gendm  `GenVisTau_decayMode` — decay mode of each generator-level visible tau.
/// @return Vector of integer decay-mode/genmatch codes per reco tau
ROOT::RVecF reco2gm(const ROOT::RVecF& gm, const ROOT::RVecF& genidx, const ROOT::RVecF& gendm) {
  ROOT::RVecI newgendm(gm.size()); // dimension = Tau (reco)
  for (size_t i=0; i<gm.size(); ++i) { // loop over reco taus
    if (gm[i]==1 or gm[i]==3) // prompt or tauonic electron
      newgendm[i] = 16;
    else if (gm[i]==2 or gm[i]==4) // prompt or tauonic muon
      newgendm[i] = 17; 
    else if (gm[i]==0) // q/g jet
      newgendm[i] = 18;
    else if (gm[i]==5) { // real generated tau => split into decay modes
      if (0<=genidx[i] and genidx[i]<gendm.size()) // index in range
        newgendm[i] = gendm[genidx[i]];
      else // presumably no match between reco Tau and GenVisTau
        newgendm[i] = gendm[genidx[i]];
    }
    else // default: unrecognized genmatch (should never happen, default gm is 0)
      newgendm[i] = -2;
  }
  return newgendm;
}

/// @brief Map each generator-level visible tau to its reconstructed decay mode.
/// @details Iterates over reconstructed taus and, for each one that is matched
///   to a GenVisTau (valid index in [0, ngen)), stores the reconstructed/regressed
///   decay mode at the corresponding generator-level position. Unmatched generator
///   taus keep the default value -1.
/// @param ngen   Number of GenVisTau objects (`nGenVisTau`).
/// @param genidx `Tau_genPartIdx` — index into the GenVisTau collection for each reco tau.
/// @param recodm `Tau_decayMode`  — reconstructed (or regressed) decay mode per reco tau.
/// @return Vector of length `ngen` with the reconstructed decay mode for each GenVisTau,
///         or -1 where no matching reco tau was found.
ROOT::RVecF gen2recodm(const float ngen, const ROOT::RVecI& genidx, const ROOT::RVecI& recodm) {
  ROOT::RVecI newrecodm(ngen,-1); // dimension = GenVisTau, default (no reco match) = -1
  for (size_t ir=0; ir<genidx.size(); ++ir) { // loop over reconstructed taus
    if (0<=genidx[ir] and genidx[ir]<ngen) {
      newrecodm[genidx[ir]] = recodm[ir];
    }
  }
  return newrecodm;
}

/// @brief Compute the reconstructed-to-generator pT ratio for each GenVisTau.
/// @details For every reconstructed tau with a valid match to a generator-level
///   visible tau the ratio ptreco / ptgen is stored at the generator-level index position.
///   Entries without a reco match, or where the generator pT is zero, keep the default value -1.
/// @param genpt  `GenVisTau_pt`   — generator-level visible-tau transverse momenta.
/// @param recopt `Tau_pt`         — reconstructed tau transverse momenta.
/// @param genidx `Tau_genPartIdx` — index into the GenVisTau collection for each reco tau.
/// @return Vector of length `nGenVisTau` containing the pT ratio for each generator tau,
///         or -1 where no valid reco match exists or the generator pT is zero.
ROOT::RVecF ptres(const ROOT::RVecF& genpt, const ROOT::RVecF& recopt, const ROOT::RVecI& genidx) {
  ROOT::RVecF ptratio(genpt.size(),-1.f); // dimension = GenVisTau, default (no reco match) = -1
  for (size_t ir=0; ir<genidx.size(); ++ir) { // loop over reconstructed taus
    const int ig = genidx[ir];
    if (ig < 0) continue;
    const size_t uig = static_cast<size_t>(ig); // convert signed int to unsigned size_t
    if (uig<genpt.size() and genpt[uig]!=0) {
      ptratio[uig] = recopt[uig]/genpt[uig];
    }
  }
  return ptratio;
}

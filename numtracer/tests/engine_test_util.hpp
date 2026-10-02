/// @file engine_test_util.hpp
/// @brief Shared helpers for the engine tests, which drive the contraction the way generated code
///        does: integer labels, explicit component tables and explicit projector denominators.
///
/// User code builds the same objects through @ref numtracer::Frame and the typed builders; these
/// helpers exist so a test can state an exact engine input (a chosen atom id, a hand-made
/// denominator) without going through that layer.
#pragma once

#include "numtracer/numtracer.hpp"

#include <string>
#include <utility>
#include <vector>

namespace ntest
{
  using namespace numtracer;

  /// Symbol names `x0 … x{n-1}` for a frame of @p n anonymous symbols.
  inline std::vector<std::string> names(int n)
  {
    std::vector<std::string> r;
    for (int i = 0; i < n; ++i) r.push_back("x" + std::to_string(i));
    return r;
  }

  // ── single Lorentz factors on integer labels (a LorentzTerm's element list) ──────────────────────
  inline LorentzFactor fmet(int a, int b) { return metric(LorentzIndex{a}, LorentzIndex{b}).front().e.front(); }
  inline LorentzFactor fvec(int a, Vlc k) { return vec(LorentzIndex{a}, Momentum{std::move(k)}).front().e.front(); }
  /// A projector factor with explicit atom ids (the user builders leave them to the Frame).
  inline LorentzFactor fproj(LorentzFactor::Kind kind, int a, int b, Vlc k, int atom, int atomS)
  {
    LorentzFactor f = projT(LorentzIndex{a}, LorentzIndex{b}, Momentum{std::move(k)}).front().e.front();
    f.kind = kind;
    f.atom = atom;
    f.atomS = atomS;
    return f;
  }
  inline LorentzFactor fprojT(int a, int b, Vlc k, int atom) { return fproj(LorentzFactor::ProjT, a, b, std::move(k), atom, -1); }
  inline LorentzFactor fprojL(int a, int b, Vlc k, int atom) { return fproj(LorentzFactor::ProjL, a, b, std::move(k), atom, -1); }
  inline LorentzFactor fprojE(int a, int b, Vlc k, int atom, int atomS)
  {
    return fproj(LorentzFactor::ProjE, a, b, std::move(k), atom, atomS);
  }
  inline LorentzFactor fprojM(int a, int b, Vlc k, int atomS) { return fproj(LorentzFactor::ProjM, a, b, std::move(k), -1, atomS); }
  inline LorentzFactor feps(int a, int b, int c, int d)
  {
    return epsilon(LorentzIndex{a}, LorentzIndex{b}, LorentzIndex{c}, LorentzIndex{d}).front().e.front();
  }

  // ── the engine entry points with explicit tables ──────────────────────────────────────────────────
  inline Poly contract(const Frame &F, const DiracChain &chain, const LorentzNet &lor,
                       const std::vector<std::array<Poly, 4>> &comp, const std::vector<Poly> &atomDen)
  {
    return ndetail::contract(F.nsym(), chain, lor, comp, atomDen, F.units());
  }
  inline DPoly contract_dressed(const Frame &F, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                const LorentzNet &lor, const std::vector<std::array<Poly, 4>> &comp,
                                const std::vector<Poly> &atomDen)
  {
    return ndetail::contract_dressed(F.nsym(), chain, slots, lor, comp, atomDen, F.units());
  }
  inline Poly contract_structural(const Frame &F, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                  const LorentzNet &lor, const std::vector<std::array<Poly, 4>> &comp,
                                  const std::vector<Poly> &atomDen)
  {
    return ndetail::contract_structural(F.nsym(), chain, slots, lor, comp, atomDen, F.units());
  }
  inline std::vector<Poly> collect_atom_denoms(const Frame &F, const std::vector<LorentzNet> &lors,
                                               const std::vector<std::array<Poly, 4>> &comp)
  {
    return ndetail::collect_atom_denoms(F.nsym(), lors, comp);
  }

  // ── SU(N) factors on integer labels (group 0, as generated code builds them) ─────────────────────
  inline SUNFac sunF(int g, int a, int b, int c) { return {SUNFacKind::F, g, a, b, c, {}}; }
  inline SUNFac sunT(int g, int a, int i, int j) { return {SUNFacKind::T, g, a, i, j, {}}; }
  inline SUNFac sunDeltaAdj(int g, int a, int b) { return {SUNFacKind::DeltaAdj, g, a, b, -1, {}}; }
  inline SUNFac sunDeltaFund(int g, int i, int j) { return {SUNFacKind::DeltaFund, g, i, j, -1, {}}; }
  inline SUNFac sunDiagFund(int g, int i, int j, std::vector<int> c2d)
  {
    return {SUNFacKind::DiagFund, g, i, j, -1, std::move(c2d)};
  }
  inline SUNFac sunDiagAdj(int g, int a, int b, std::vector<int> c2d)
  {
    return {SUNFacKind::DiagAdj, g, a, b, -1, std::move(c2d)};
  }

  // ── whole-net builders on integer labels and a single momentum id, as generated code writes them ──
  inline LorentzNet imet(int a, int b) { return leaf({.kind = LorentzFactor::Metric, .a = a, .b = b}); }
  inline LorentzNet ivec(int a, int vid) { return leaf({.kind = LorentzFactor::Vector, .a = a, .b = -1, .vlc = {{1.0, vid}}}); }
  inline LorentzNet iprojT(int a, int b, int vid, int atom)
  {
    return leaf({.kind = LorentzFactor::ProjT, .a = a, .b = b, .vid = vid, .atom = atom});
  }
  inline LorentzNet iprojL(int a, int b, int vid, int atom)
  {
    return leaf({.kind = LorentzFactor::ProjL, .a = a, .b = b, .vid = vid, .atom = atom});
  }
  inline LorentzNet iprojE(int a, int b, int vid, int atom, int atomS)
  {
    return leaf({.kind = LorentzFactor::ProjE, .a = a, .b = b, .vid = vid, .atom = atom, .atomS = atomS});
  }
  inline LorentzNet iprojM(int a, int b, int vid, int atomS)
  {
    return leaf({.kind = LorentzFactor::ProjM, .a = a, .b = b, .vid = vid, .atomS = atomS});
  }
  inline LorentzNet iepsilon(int a, int b, int c, int d)
  {
    return leaf({.kind = LorentzFactor::Epsilon, .a = a, .b = b, .c = c, .d = d});
  }
} // namespace ntest

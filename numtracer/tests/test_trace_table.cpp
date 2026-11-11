/// @file test_trace_table.cpp
/// @brief Trace tables (codegen/trace_table.hpp): the table form of real projector traces reproduces the
///        traces, the lowered fill program reproduces the coefficients, dressed and complex traces are
///        tabulated piece by piece, the balance keeps only the tables that pay, and wrong scale symbols are caught.
///
/// Frame: an external scale S with a shape angle (p1 = S e1, p2 = S (cos chi, sin chi, 0, 0)) and a loop
/// momentum l = L (c, s cos phi, s sin phi, 0). The table axes are the angles; L and S are the scales.

#include "numtracer/numtracer.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace nt = numtracer;

namespace
{
  int fails = 0;
  void check(bool ok, const std::string &what)
  {
    std::printf("  %-66s %s\n", what.c_str(), ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }

  std::string sci(double v)
  {
    char b[32];
    std::snprintf(b, sizeof b, "%.1e", v);
    return b;
  }

  bool throws(const std::function<void()> &f)
  {
    try {
      f();
    } catch (const std::exception &) {
      return true;
    }
    return false;
  }
} // namespace

int main()
{
  std::mt19937 rng(20261111);
  std::uniform_real_distribution<double> C(-0.95, 0.95), U(0.0, 1.0);

  nt::Frame F;
  auto S = F.symbol("S"), L = F.symbol("L");
  auto [Cc, Sc] = F.angle("chi");
  auto [Ct, St] = F.angle("theta");
  auto [Cp, Sp] = F.angle("phi");
  auto p1 = F.momentum(S, 0, 0, 0);
  auto p2 = F.momentum(S * Cc, S * Sc, 0, 0);
  auto l = F.momentum(L * Ct, L * St * Cp, L * St * Sp, 0);
  auto [mu, nu, a, b] = F.indices<4>();

  struct Case {
    std::string name;
    nt::Poly T;
    nt::TableScales sc;
  };
  std::vector<Case> cases;
  cases.push_back({"tr[p1/ g^mu l/ g^nu] P_T(l-p1)",
                   F.trace({nt::slash(p1), nt::gamma(mu), nt::slash(l), nt::gamma(nu)}, nt::projT(mu, nu, l - p1)),
                   {L.id, S.id}});
  cases.push_back({"tr[p1/ g l/ g p2/ g l/ g] P_T(l-p1) P_T(l+p2)",
                   F.trace({nt::slash(p1), nt::gamma(mu), nt::slash(l), nt::gamma(nu), nt::slash(p2), nt::gamma(a),
                            nt::slash(l), nt::gamma(b)},
                           nt::projT(mu, a, l - p1) * nt::projT(nu, b, l + p2)),
                   {L.id, S.id}});
  cases.push_back({"p2.P_T(l-p2).p1",
                   F.contract(nt::vec(mu, p2) * nt::projT(mu, nu, l - p2) * nt::vec(nu, p1)), {L.id, S.id}});
  cases.push_back({"l.l (no external scale)", F.contract(nt::vec(mu, l) * nt::vec(mu, l)), {L.id, -1}});

  auto point = [&](double s, double lv) {
    return F.at({{S, s}, {L, lv}, {Cc, C(rng)}, {Ct, C(rng)}, {Cp, C(rng)}});
  };
  auto inv_values = [&](const nt::Point &pt) {
    std::vector<double> v(F.denominators().size());
    for (std::size_t i = 0; i < v.size(); ++i)
      v[i] = F.eval(F.inv_denominator(static_cast<int>(i)), pt).re;
    return v;
  };

  std::printf("== table form reproduces the trace ==\n");
  std::vector<nt::TableTrace> tables;
  for (const Case &c : cases) {
    const auto t = nt::table_form(c.T, F.denominators(), c.sc);
    check(t.has_value(), c.name + ": tabulable");
    if (!t) continue;
    std::printf("    d=%d D=%d n=%d..%d entries=%d atoms=%zu\n", t->d, t->D, t->n0, t->n1, t->entries(), t->inv.size());
    // l/S across 2 decades: the direct evaluation is still accurate there, so the comparison is tight
    double worst = 0;
    for (int k = 0; k < 200; ++k) {
      const double s = 0.3 + 2 * U(rng), lv = s * std::pow(10.0, -1 + 2 * U(rng));
      const nt::Point pt = point(s, lv);
      std::vector<double> av;
      for (const nt::Poly &ac : t->coeff)
        av.push_back(F.eval(ac, pt).re);
      const double ref = F.eval(c.T, pt).re;
      const double got = nt::eval_table_trace(*t, av.data(), lv, s, inv_values(pt));
      worst = std::max(worst, std::fabs(got - ref) / std::max(1.0, std::fabs(ref)));
    }
    check(worst < 1e-10, c.name + ": reconstruction = trace (worst " + sci(worst) + ")");
    tables.push_back(*t);
  }

  std::printf("== the lowered fill reproduces the coefficients ==\n");
  {
    std::vector<const nt::TableTrace *> ptr;
    for (const auto &t : tables)
      ptr.push_back(&t);
    const nt::TableFillProg prog = nt::lower_table_fill(ptr);
    std::size_t nEntries = 0;
    for (const auto &t : tables)
      nEntries += static_cast<std::size_t>(t.entries());
    check(prog.roots.size() == nEntries, "one root per entry");
    double worst = 0;
    for (int k = 0; k < 50; ++k) {
      const nt::Point pt = point(1.0, 1.0);
      const std::vector<double> f = F.fill_values(prog.env, pt);
      std::size_t e = 0;
      for (const auto &t : tables)
        for (const nt::Poly &ac : t.coeff) {
          const double ref = F.eval(ac, pt).re;
          const double got = nt::interpret(prog.ins, prog.roots[e++], f.data());
          worst = std::max(worst, std::fabs(got - ref) / std::max(1.0, std::fabs(ref)));
        }
    }
    check(worst < 1e-12, "interpret(fill) = coefficients (worst " + sci(worst) + ")");

    std::ostringstream os;
    nt::FillFormulas fm;
    fm.var = [&](int id) { return F.symbol_names().at(static_cast<std::size_t>(id)); };
    nt::emit_table_fill(os, prog, fm, "tab_fill", "double cos_chi, double cos_theta, double cos_phi", "inline");
    check(os.str().find("void tab_fill(double *out, const std::size_t stride,") != std::string::npos &&
              os.str().find("out[" + std::to_string(nEntries - 1) + " * stride]") != std::string::npos,
          "emit_table_fill: signature and every entry written");
    check(os.str().find("float") == std::string::npos, "emit_table_fill: fp64 throughout");
  }

  std::printf("== dressed and complex traces ==\n");
  // T = D0 (0.3 + 0.7 i) T1 + D0 D1 T2 + T3 i, with T1..T3 the projector traces above
  nt::DressedPoly dp = nt::DressedPolyFactory::zero(cases[0].T.nsym);
  dp.add({0}, cases[0].T * F.constant(nt::Cx{0.3, 0.7}));
  dp.add({0, 1}, cases[1].T);
  dp.add({}, cases[2].T * F.constant(nt::Cx{0, 1}));
  nt::GlobalEnv genv;
  const nt::TableScales sc{L.id, S.id};
  const auto dt = nt::tabulate(dp, F.denominators(), sc, /*complex=*/true, genv);
  check(dt.has_value() && dt->complex, "dressed complex trace: tabulable");
  if (dt) {
    std::size_t nIm = 0;
    for (const auto &p : dt->pieces)
      nIm += p.imag;
    check(dt->pieces.size() == 4 && nIm == 2, "one piece per (dressing monomial, re/im) part (4, 2 imaginary)");
    const std::vector<double> dress = {1.7, -0.4};
    double worst = 0;
    for (int k = 0; k < 200; ++k) {
      const double s = 0.3 + 2 * U(rng), lv = s * std::pow(10.0, -1 + 2 * U(rng));
      const nt::Point pt = point(s, lv);
      std::vector<double> av;
      for (const auto &p : dt->pieces)
        for (const nt::Poly &ac : p.t.coeff)
          av.push_back(F.eval(ac, pt).re);
      const nt::Cx ref = F.eval(dp, pt, dress);
      const nt::Cx got = nt::eval_tabulated_trace(*dt, av.data(), lv, s, inv_values(pt), dress);
      worst = std::max(worst, std::hypot(got.re - ref.re, got.im - ref.im) / std::max(1.0, std::hypot(ref.re, ref.im)));
    }
    check(worst < 1e-10, "reconstruction = dressed complex trace (worst " + sci(worst) + ")");
  }
  const auto rt = nt::tabulate(dp, F.denominators(), sc, /*complex=*/false, genv);
  check(rt && !rt->complex && rt->pieces.size() == 2, "a real-consumed trace keeps only its real parts");
  {
    nt::GlobalEnv cg;
    cg.cplxDress = {0, 1};
    check(!nt::tabulate(dp, F.denominators(), sc, true, cg), "a complex dressing -> not tabulable");
  }

  std::printf("== the balance ==\n");
  {
    std::vector<std::optional<nt::TabulatedTrace>> tabs = {dt, dt, std::nullopt};
    std::vector<nt::GenProg> progs(3);
    progs[0].ins.resize(static_cast<std::size_t>(dt->recon_ops() + 100 * dt->entries())); // 100 instr per entry
    progs[1].ins.resize(static_cast<std::size_t>(dt->recon_ops() + 2 * dt->entries()));   // 2 per entry
    nt::GlobalEnv g;
    auto sel = tabs;
    nt::select_table_traces(sel, progs, 1.0, sc, g); // fp64 emission: 8 bytes per entry
    check(sel[0] && !sel[1] && !sel[2], "balance 1: keeps 100 instr/entry, drops 2 instr/entry");
    sel = tabs;
    nt::select_table_traces(sel, progs, 0.0, sc, g);
    check(sel[0] && sel[1], "balance 0: keeps every table");
  }

  std::printf("== header emission ==\n");
  {
    nt::GlobalEnv g;
    std::vector<std::optional<nt::TabulatedTrace>> tabs = {dt, std::nullopt, dt};
    std::vector<nt::GenProg> progs(3);
    nt::select_table_traces(tabs, progs, 0.0, sc, g);
    nt::FillFormulas fm;
    fm.var = [&](int id) { return F.symbol_names().at(static_cast<std::size_t>(id)); };
    std::ostringstream os;
    nt::emit_trace_tables(os, tabs, {0, 1, 0}, fm, "double cos_chi, double cos_theta, double cos_phi",
                          {"cos_chi", "cos_theta", "cos_phi"}, sc, g, "static inline");
    const std::string s = os.str();
    auto has = [&](const std::string &x) { return s.find(x) != std::string::npos; };
    check(has("void tab_fill0(double *out, const std::size_t stride,") &&
              has("tab_fill0(out + 0 * stride, stride, cos_chi, cos_theta, cos_phi);") && !has("tab_fill2("),
          "one fill per distinct tabulated trace, called from tab_fill");
    check(has("nt_complex_t ttr0(") && has("return nt_complex_t{re_, im_};") && has("fma(h_, r_,") &&
              has("h_ *= f[") && has("h_ = tab_scale<"),
          "ttr0: complex, Horner in r, dressing factors, s^d applied last");
    check(has("ttr1([[maybe_unused]] const double *row, [[maybe_unused]] const std::size_t st, const double *f) { return tr1(f); }"),
          "ttr1: inline forwarder");
    check(has("ttr2(const double *row, const std::size_t st, const double *f) { return ttr0(row, st, f); }"),
          "ttr2: duplicate of trace 0, forwarded");
    check(has("tab_size = " + std::to_string(dt->entries()) + ";"), "tab_size counts the distinct table once");
  }

  std::printf("== guards ==\n");
  check(throws([&] { (void)nt::table_form(cases[0].T, F.denominators(), {L.id, -1}); }),
        "missing external scale -> not homogeneous -> throws");
  check(throws([&] { (void)nt::table_form(cases[0].T, F.denominators(), {Ct.id, S.id}); }),
        "an angle as the radial scale -> throws");
  check(throws([&] { (void)nt::table_form(cases[0].T, F.denominators(), {-1, S.id}); }), "no radial scale -> throws");

  std::printf(fails ? "\n%d FAILED\n" : "\nall passed\n", fails);
  return fails ? 1 : 0;
}

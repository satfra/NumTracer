// step-23 — Trace tables: fill once per run, evaluate at every scale.
//
// trace_tables.wls emitted the same network twice: Tab_kernel with every trace tabulated, and
// Tab_auto_kernel with NumTracer deciding per trace. This program does what an integrator does with
// a table: it fills one row per angle node, once, and then evaluates kernel_tab() at many (l1, p).
// Both entry points must agree with the closed form; the table is the more accurate of the two.
#include "Tab_auto_kernel.hh"
#include "Tab_kernel.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
int failures = 0;

// Any callable works as a dressing; the table never sees it.
struct MyZ {
  double operator()(double x) const { return 1.0 / (1.0 + x * x); }
};

// The closed form of the network (see trace_tables.wls).
double closed(double l, double c, double p, const MyZ &Z) {
  const double transverse = l * p * p * (1 - c * c) * (l - c * p) / (l * l - 2 * c * l * p + p * p);
  const double dirac = 4 * l * l * p * p * (2 * c * c - 1) - 4 * l * p * p * p * c;
  return transverse + Z(l) * dirac;
}

void check(const char *what, bool ok) {
  if (!ok)
    ++failures;
  std::printf("  %-52s %s\n", what, ok ? "ok" : "FAILED");
}
} // namespace

int main() {
  const MyZ myZ;

  // @snip begin: traits
  // What the kernel class declares about its table. tab_box_mask has one entry per integration
  // variable (l1, cos1) and tab_ext_mask one per external coordinate (p): true marks a variable the
  // table depends on, i.e. a table axis. The scales l1 and p are not axes.
  static_assert(Tab_kernel::has_tab);
  std::printf("Tab_kernel: tab_size = %zu entries per node, axes: l1 %d, cos1 %d, p %d\n", Tab_kernel::tab_size,
              Tab_kernel::tab_box_mask[0], Tab_kernel::tab_box_mask[1], Tab_kernel::tab_ext_mask[0]);
  // @snip end: traits

  // @snip begin: fill
  // The table: tab_size entries for each of n nodes, node-fastest. Entry e of node i sits at
  // table[e * n + i], so node i's row is table.data() + i with stride n. tab() takes kernel()'s
  // arguments and reads only the axes, so the scales and the dressing may hold any value.
  const std::size_t n = 10;
  std::vector<double> nodes(n), table(Tab_kernel::tab_size * n);
  for (std::size_t i = 0; i < n; ++i) {
    nodes[i] = -1.0 + (2.0 * i + 1.0) / n; // cos1 = -0.9, -0.7, ..., 0.9
    Tab_kernel::tab(table.data() + i, n, /*l1*/ 1.0, nodes[i], /*p*/ 1.0, myZ);
  }
  // @snip end: fill

  std::printf("the entries of the node cos1 = %g:", nodes[7]);
  for (std::size_t e = 0; e < Tab_kernel::tab_size; ++e)
    std::printf(" %g", table[e * n + 7]);
  std::printf("\n");

  // The same fill with other scales must give the same table, bit for bit.
  std::vector<double> other(table.size());
  for (std::size_t i = 0; i < n; ++i)
    Tab_kernel::tab(other.data() + i, n, /*l1*/ 37.0, nodes[i], /*p*/ 0.01, myZ);
  const bool sameFill = other == table;

  // @snip begin: use
  // Evaluation: kernel_tab() takes the node's row and the stride in front of kernel()'s arguments.
  // One table serves every l1 and p, on both sides of l1 = p. Both entry points are compared with
  // the closed form, relative to the largest |value| over the nodes.
  std::printf("relative deviation from the closed form, worst over the %zu nodes\n", n);
  std::printf("  %6s %8s %12s %12s\n", "p", "l1/p", "kernel_tab", "kernel");
  double worstTab = 0, worstInline = 0;
  for (const double p : {0.1, 1.0, 30.0})
    for (const double ratio : {1e-3, 0.1, 1.0, 10.0, 1e3}) {
      const double l1 = ratio * p;
      double scale = 0, dTab = 0, dInline = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const double ref = closed(l1, nodes[i], p, myZ);
        const double tab = Tab_kernel::kernel_tab(table.data() + i, n, l1, nodes[i], p, myZ);
        const double inl = Tab_kernel::kernel(l1, nodes[i], p, myZ);
        scale = std::max(scale, std::fabs(ref));
        dTab = std::max(dTab, std::fabs(tab - ref));
        dInline = std::max(dInline, std::fabs(inl - ref));
      }
      std::printf("  %6g %8g %12.1e %12.1e\n", p, ratio, dTab / scale, dInline / scale);
      worstTab = std::max(worstTab, dTab / scale);
      worstInline = std::max(worstInline, dInline / scale);
    }
  // @snip end: use

  std::printf("checks\n");
  check("the fill does not depend on l1 or p", sameFill);
  check("kernel_tab == closed form (1e-14)", worstTab < 1e-14);
  check("kernel == closed form (1e-11)", worstInline < 1e-11);

  // @snip begin: auto
  // Under Automatic both traces are cheaper inline than read from a table, so nothing is tabulated:
  // tab_size is 0, has_tab is false and kernel_tab() ignores its row. An integrator then simply
  // calls kernel().
  static_assert(!Tab_auto_kernel::has_tab && Tab_auto_kernel::tab_size == 0);
  const double viaTab = Tab_auto_kernel::kernel_tab(nullptr, 0, 0.7, 0.3, 1.9, myZ);
  const double inl = Tab_auto_kernel::kernel(0.7, 0.3, 1.9, myZ);
  check("Automatic: nothing tabulated, kernel_tab == kernel", std::fabs(viaTab - inl) <= 1e-14 * std::fabs(inl));
  // @snip end: auto

  std::printf(failures == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return failures == 0 ? 0 : 1;
}

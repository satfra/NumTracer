# 06 — Generating kernels from Mathematica

This chapter has no compiled program: it is a walkthrough, not a hand-coded exercise. Tutorials
1–5 built the quark self-energy by hand in C++; the Mathematica front-end automates exactly that
— analyse the network, fix the frame, contract each diagram numerically, lower it — and writes a
flat C++ kernel.

See the narrative page `numtracer/documentation/tutorials/generating-kernels.md`. The driver is
`numtracer/tests/gen/gen_zq_numeric.wls` (`NumTrace` + `MakeNTKernel "Backend" -> "Numeric"`), and
the committed artifact it produces is `numtracer/tests/gen/Zq_num_kernel.hh`.

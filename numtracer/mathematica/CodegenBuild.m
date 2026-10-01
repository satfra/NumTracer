(* per-job RAM estimate (GB) used to bound the parallel compile. With defs chunked, a unit's peak RSS
   is a few hundred MB, so this is deliberately conservative headroom. NT_GEN_JOBS pins the job
   count outright. *)

$ntGenJobMemGB = 1.5;

(* MemAvailable: what the kernel can hand out without swapping. Read fresh (not at package load) —
   the Wolfram kernel that just built the flow may itself be holding many GB. *)

ntAvailMemGB[] := Quiet @
    Check[
      Module[{m = StringCases[Import["/proc/meminfo", "Text"], "MemAvailable:" ~~ Whitespace ~~ d : DigitCharacter.. ~~ Whitespace ~~ "kB" :> ToExpression[d], 1]},
        If[m === {},
          N[MemoryAvailable[] / 2^30],
          First[m] / 1024. / 1024.]],
      $Failed];

(* max concurrent cxx processes in the generator compile. Peak RAM ~ jobs x per-unit RSS, so this is
   bounded by BOTH cores and free memory: a cores-only cap will OOM the box on a dense flow, where a
   single unit can take GBs. Def chunking now bounds a unit's size, so the memory term is belt-and-braces
   that keeps a future monster flow from thrashing instead of failing loudly. NT_GEN_JOBS pins the count
   outright. Delayed (:=) so the memory reading is taken at compile time, not at package load. *)

ntCompileJobs[] := With[{e = ntEnvPosInt["NT_GEN_JOBS"], avail = ntAvailMemGB[]},
    Which[
      e > 0,
        e,
      IntegerQ[$ProcessorCount],
        Max[
          2,
          Min[
            24,
            $ProcessorCount,
            If[NumericQ[avail] && avail > 0,
              Floor[avail / $ntGenJobMemGB],
              24]]],
      True,
        4]];

$ntCompileJobs := ntCompileJobs[];

(* C++ compiler for the build-time generator. The emitted generator is ordinary numeric C++, so
   either g++ or clang++ compiles it with the same flags (-std=c++20 -ftemplate-depth=4000
   -O1/-O0 -fno-exceptions -fno-rtti -pthread -I -c -o). Resolution order: the NT_GEN_CXX env
   override (verbatim), else PREFER clang++ when present (it compiles the -O0 net-builder units
   ~2.5x faster than g++ — u33 5.8 s -> 2.3 s), else a compiler detected via CCompilerDriver
   (GCC -> g++, Clang -> clang++; binary taken from the driver's CompilerInstallation dir), else
   "g++" on PATH. clang objects link fine against the g++-built libNumTracer.a (shared libstdc++
   ABI on Linux). *)

resolveGenCxx[] := Module[{env = Environment["NT_GEN_CXX"], path, comps, clangPick, pick, name, inst, exe, full},
    If[StringQ[env] && StringTrim[env] =!= "",
      Return[StringTrim[env]]];
    Quiet @ Needs["CCompilerDriver`"];
    comps = Quiet @ Check[CCompilers[], {}];
    If[!ListQ[comps],
      comps = {}];
    comps = Select[comps, AssociationQ[#] && StringQ[Lookup[#, "CompilerInstallation"]]&];
    (* 1) prefer a Clang known to CCompilerDriver *)
    clangPick = SelectFirst[comps, StringContainsQ[ToString @ Lookup[#, "Name", ""], "Clang", IgnoreCase -> True]&, None];
    If[clangPick =!= None,
      full = FileNameJoin[{Lookup[clangPick, "CompilerInstallation"], "clang++"}];
      If[FileExistsQ[full], Return[full]]];   (* else: clang++ on PATH, or the g++ fallback below *)
    (* 2) prefer clang++ on PATH even if CCompilerDriver did not enumerate it *)
    path = Environment["PATH"];
    If[StringQ[path] && AnyTrue[StringSplit[path, ":"], FileExistsQ[FileNameJoin[{#, "clang++"}]]&],
      Return["clang++"]];
    (* 3) fall back to the CCompilerDriver default (typically g++) *)
    If[comps === {},
      Return["g++"]];
    pick = SelectFirst[comps, Lookup[#, "Compiler"] === Quiet @ DefaultCCompiler[]&, First[comps]];
    name = ToString @ Lookup[pick, "Name", ""];
    inst = Lookup[pick, "CompilerInstallation"];
    exe =
      If[StringContainsQ[name, "Clang", IgnoreCase -> True],
        "clang++",
        "g++"];
    full = FileNameJoin[{inst, exe}];
    If[FileExistsQ[full],
      full,
      exe]];

(* tcmalloc preload for the generator RUN: the trace engine's phase A/B does heavy tiny-alloc
   churn and a tcmalloc_minimal LD_PRELOAD is a measured -7% on the run (it also lowers the
   glibc-arena fragmentation floor, see PHASE-A residue notes). Env prefix on the Run[] shell
   command, so it needs no code in the generator itself. Silently empty when the library is absent,
   so nothing changes on machines without gperftools — which is also the opt-out. *)
ntTcmallocPrefix[] :=
  With[{lib = SelectFirst[
      {"/usr/lib/libtcmalloc_minimal.so", "/usr/lib64/libtcmalloc_minimal.so",
       "/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4", "/usr/lib/libtcmalloc_minimal.so.4"},
      FileExistsQ, None]},
    If[lib === None, "", "LD_PRELOAD='" <> lib <> "' "]];

(* NT_GEN_DEVICE=1 tells the generator that the emitted kernel is DEVICE code, which is what enables
   gen.hpp's size-gated `__noinline__` (device-only: the host has no register cliff, and its
   all-inline emission is byte-identical). It is passed explicitly rather than sniffed from the
   decorator because `ntKokkosDecor` rewrites the raw CUDA spelling to the Kokkos macros before this
   point — that is exactly why the old `decor.find("__device__")` test silently never fired on any
   production flow. `Automatic` falls back to the raw spelling so existing raw-CUDA callers keep
   working; DiFfRG_compat passes True/False from its own Device option. *)
(* The single definition of "does this emission target device code". Used twice: for the ONLINE
   Run[] env prefix below, and for the manifest's "device" field, which is how the OFFLINE numtrace
   build step learns the same fact (it inherits nothing of this kernel's environment — see the
   identical argument for the thread caps in ntWriteManifest). Keeping one predicate is the point:
   the two paths must never disagree about whether the gate applies. *)
ntDeviceTargetQ[deviceTarget_, decor_String] :=
  If[deviceTarget === Automatic,
    StringContainsQ[decor, "__device__"],
    TrueQ[deviceTarget]];

ntDeviceEnvPrefix[deviceTarget_, decor_String] :=
  If[ntDeviceTargetQ[deviceTarget, decor], "NT_GEN_DEVICE=1 ", ""];

(* ---- offline generation: manifest + verdict plumbing --------------------------------------
   The verdict macro for a flow's C++ namespace tag: za_qcd -> NT_ZA_QCD_VERDICT. *)
ntVerdictMacro[ns_String] := "NT_" <> ToUpperCase[StringReplace[ns, Except[WordCharacter] -> "_"]] <> "_VERDICT";
ntVerdictFile = "numtrace_verdict.hh";
ntManifestFile = "numtrace.json";

(* Canonicalise: resolve symlinks so two spellings of the SAME directory compare equal. Without this,
   a project reached through a symlinked home (/home/me/Code -> /mnt/data/Code) yields a flow dir and a
   gen dir that share no prefix, and the "relative" path becomes a deep ../../.. climb to the root and
   back down — resolvable but machine-specific and unusable in a committed manifest.
   AbsoluteFileName resolves links but requires the path to EXIST. Offline, the traces header does not
   exist yet when the manifest is written, so a naive fallback to ExpandFileName left it spelled
   /home/... while its (existing) parent directory resolved to /mnt/data/... — two spellings of the
   same place with no common prefix, which made the "relative" path degenerate to an absolute one and
   the build then joined it onto the flow dir again (…/flows/ZA3//home/…/flows/ZA3/kernels.hh). So
   canonicalise the deepest ANCESTOR that does exist and re-append the segments below it. *)
ntCanonicalDir[p_String] := Module[{a = Quiet@Check[AbsoluteFileName[p], $Failed], parent, base},
  If[StringQ[a], Return[a]];
  base = FileNameTake[p];
  parent = FileNameDrop[p, -1];
  (* no progress possible (root, or a bare relative name): fall back to plain expansion *)
  If[parent === "" || parent === p || base === "", Return[ExpandFileName[p]]];
  FileNameJoin[{ntCanonicalDir[parent], base}]];

ntRelativePath[from_String, to_String] := Module[{f, t, common, rel},
  f = DeleteCases[FileNameSplit[ntCanonicalDir[from]], ""];
  t = DeleteCases[FileNameSplit[ntCanonicalDir[to]], ""];
  common = LengthWhile[Range[Min[Length[f], Length[t]]], f[[#]] === t[[#]] &];
  (* A path that still has to climb out to the root shares nothing meaningful with the target
     (different mount, say). An absolute path is at least honest about that. *)
  If[common <= 1 && Length[f] > 1, Return[ntCanonicalDir[to]]];
  rel = FileNameJoin[Join[ConstantArray["..", Length[f] - common], Drop[t, common]]];
  If[rel === "", ".", rel]];

(* ---- per-flow numtrace manifest -------------------------------------------------------------
   Written beside the kernel headers (mirroring DiFfRG's per-flow sources.m, aggregated the same way)
   and COMMITTED with them. Two jobs: the "generated" 0|1 switch — emitting sets 0, the numtrace target
   sets 1 once the kernels are actually built (a fresh clone has the committed kernels but no gen/, so
   timestamps cannot express "already generated") — and the metadata the build needs (namespace, source
   list, decorator, main-TU -O level, whether a probe is required). One writer per file, so parallel
   numtrace jobs never race. Paths are basenames relative to "gen_dir", relative to the flow dir. *)


ntWriteManifest[flowDir_String, spec_Association] := Module[
  {name = spec["Class"], ns = spec["Namespace"], genFile = spec["Generator"], tracesFile = spec["Traces"],
   unitFiles = spec["Units"], decor = spec["Decorator"], mainOpt = spec["MainOpt"],
   complexQ = spec["Complex"], probeFile = spec["Probe"], deviceTarget = spec["DeviceTarget"],
   genDir, manifest},
  genDir = DirectoryName[genFile];
  manifest = <|
    (* the flow's identity is its directory (flows/ZA4), not the kernel class name (ZA4_kernel) *)
    "name"          -> FileNameTake[StringTrim[flowDir, "/"]],
    "class"         -> name,
    "namespace"     -> ns,
    "generated"     -> 0,
    "gen_dir"       -> ntRelativePath[flowDir, genDir],
    "generator"     -> FileNameTake[genFile],
    "units"         -> FileNameTake /@ unitFiles,
    "decorator"     -> decor,
    "main_opt"      -> mainOpt,
(* the per-flow thread caps, captured from whatever SetNumTracerThreads[nA, nB] was in force at emit
   time — it sets exactly these two environment variables. Offline the generator runs from a cmake -P
   build step, which inherits nothing of the emitting Wolfram kernel's environment, so a cap that is
   not written down here is simply lost. The numtrace driver re-applies them, lowering (never raising)
   the build's own -jN. 0 = unset. *)
    "maxw"          -> ntEnvPosInt["NT_GEN_MAXW"],
    "maxw_b"        -> ntEnvPosInt["NT_GEN_MAXW_B"],
(* Does this flow's kernel target DEVICE code? Same trip, same reason as the thread caps above: it
   enables gen.hpp's size-gated `__noinline__`, which is device-only, and the generator learns it
   from NT_GEN_DEVICE. Online that variable comes from ntDeviceEnvPrefix's shell prefix; offline
   there is no shell prefix and no inherited environment, so without this field the gate is simply
   dead — which is exactly what it was on every offline-generated flow until this was added. *)
    "device"        -> ntDeviceTargetQ[deviceTarget, decor],
    "complex"       -> TrueQ[complexQ],
    "kernels"       -> ntRelativePath[flowDir, tracesFile]|>;
  If[TrueQ[complexQ] && StringQ[probeFile],
    manifest = Join[manifest, <|
      "probe"         -> FileNameTake[probeFile],
      "verdict_macro" -> ntVerdictMacro[ns],
      "verdict"       -> ntVerdictFile|>]];
  Export[FileNameJoin[{flowDir, ntManifestFile}], manifest, "JSON"];
  FileNameJoin[{flowDir, ntManifestFile}]];

(* Flip a manifest's switch to 1 — the flow's kernels are built and committed. The offline twin of
   this lives in NumTracerNumtraceRun.cmake, run by the build once the generator has succeeded. *)
ntMarkGenerated[manifestFile_String] := Module[{m = Import[manifestFile, "RawJSON"]},
  Export[manifestFile, Append[m, "generated" -> 1], "JSON"]];

(* Locate the NumTracer C++ headers for the generator compile. In-tree the package sits in
   <repo>/numtracer/mathematica/ with the headers in ../include; installed (e.g. under
   $UserBaseDirectory/Applications/NumTracer) the header location is recorded at configure time
   in the CMake-generated sibling NumTracerPaths.m. Overridable per call ("IncludeDir" option)
   or globally (NUMTRACER_INCLUDE_DIR environment variable). *)

resolveIncludeDir::nodir = "Cannot locate the NumTracer C++ headers (tried the NUMTRACER_INCLUDE_DIR environment variable, the in-tree ../include, the installed NumTracerPaths.m record, and ~/.local/share/NumTracer/include). Pass \"IncludeDir\" -> dir to MakeNTKernel.";

resolveIncludeDir[] := Module[{envDir, dir, pathsFile},
    envDir = Environment["NUMTRACER_INCLUDE_DIR"];
    If[StringQ[envDir] && DirectoryQ[envDir],
      Return[envDir]];
    dir = FileNameJoin[{DirectoryName[$NumTracerDirectory], "include"}];
    If[DirectoryQ[dir],
      Return[dir]];
    pathsFile = FileNameJoin[{$NumTracerDirectory, "NumTracerPaths.m"}];
    If[FileExistsQ[pathsFile],
      Get[pathsFile];
      If[StringQ[$NumTracerInstalledIncludeDir] && DirectoryQ[$NumTracerInstalledIncludeDir],
        Return[$NumTracerInstalledIncludeDir]]];
    dir = FileNameJoin[{$HomeDirectory, ".local", "share", "NumTracer", "include"}];
    If[DirectoryQ[dir],
      Return[dir]];
    Message[resolveIncludeDir::nodir];
    Abort[]];

(* Locate the compiled engine library `libNumTracer.a` for the generator LINK. NumTracer ships as a
   compiled static library by default: the heavy engine bodies (numeric contraction / SU(N) fold /
   lowering) are compiled ONCE into this archive, so the emitted generator only parses the declaration
   headers and links the archive instead of re-instantiating and -O2-optimising the whole engine on
   every generation run. Searched: the NT_GEN_LIB environment variable (a full path), the installed
   `<prefix>/lib` sibling of the include dir, and the in-tree `<repo>/numtracer/build` default build
   dir. Returns the path, or $Failed — in which case the generator falls back to a self-contained
   header-only compile (`-DNUMTRACER_HEADER_ONLY=1`), which still works but pays the old compile floor. *)
(* The generator compiles gen_*.cpp against the CURRENT headers in incDir but LINKS the prebuilt
   archive. Those two must come from the same source vintage: the engine's types (MPoly's storage and
   its construction path, SUNNet, the fold buffers) are defined in the headers and compiled into the
   archive, so a header edit that is not followed by a library rebuild is a silent ODR/ABI mismatch —
   the generator RUNS, emits a plausible kernels.hh, and the numbers are wrong. Nothing downstream can
   catch it: kernel.hh only NAMES the traces (s1..sN) and is emitted by the Mathematica layer, so it
   stays byte-identical while every trace body in kernels.hh changes underneath it — which is exactly
   how a stale archive silently corrupted a flow's ZA4 and survived a "kernel.hh unchanged" check.
   mtime is a coarse proxy for provenance, but the failure it guards is silent and total, so err loud:
   a library older than any header it was built from is never trustworthy. *)

MakeNTKernel::stalelib = "The prebuilt engine archive\n  `1`\nis OLDER than the header\n  `2`\nthat the generator will compile against. Linking them mixes two source vintages of the engine's types (MPoly/SUNNet/fold buffers) — an ODR/ABI mismatch that SILENTLY produces wrong traces in kernels.hh while leaving kernel.hh byte-identical. Rebuild the library first, e.g.\n  cmake --build <repo>/numtracer/build --target NumTracer\nthen regenerate. (Set NT_GEN_LIB to a specific archive, or NT_ALLOW_STALE_LIB=1 to override — the latter is almost never right.)";

genLibStaleQ[lib_, incDir_] := Module[{hdrs, newest},
    hdrs = FileNames["*.hpp" | "*.h", incDir, Infinity];
    If[hdrs === {},
      Return[False]];
    newest = Max[AbsoluteTime /@ (FileDate[#, "Modification"]& /@ hdrs)];
    AbsoluteTime[FileDate[lib, "Modification"]] < newest];

resolveGenLib[incDir_] := Module[{env, base, cands, lib},
    env = Environment["NT_GEN_LIB"];
    lib =
      If[StringQ[env] && FileExistsQ[env],
        env,
        base = DirectoryName[incDir];(* <prefix> (installed) or <repo>/numtracer (in-tree): parent of include/ *)
        cands =
          {
            FileNameJoin[{base, "lib", "libNumTracer.a"}],
            (* installed layout *)
            FileNameJoin[{base, "build", "libNumTracer.a"}]
          };(* in-tree default build dir *)
        SelectFirst[cands, FileExistsQ, $Failed]];
    If[lib =!= $Failed && !ntEnvFlag["NT_ALLOW_STALE_LIB"] && genLibStaleQ[lib, incDir],
      Module[{hdrs = FileNames["*.hpp" | "*.h", incDir, Infinity], newest},
        newest = First @ SortBy[hdrs, -AbsoluteTime[FileDate[#, "Modification"]]&];
        Message[MakeNTKernel::stalelib, lib, newest];
        Abort[]]];
    lib];

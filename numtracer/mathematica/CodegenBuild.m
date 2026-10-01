(* CodegenBuild.m — build plumbing for the generator program: parallel-compile job count, compiler
   choice, tcmalloc/device env prefixes, include dir and engine library lookup (with the stale-archive
   guard), and the per-flow numtrace manifest. Loaded by NumTracer.m via ntLoadPart, in
   NumTracer`Private`. *)

(* per-job RAM estimate (GB) used to bound the parallel compile; deliberately conservative, since a
   chunked unit peaks at a few hundred MB. NT_GEN_JOBS pins the job count outright. *)

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

(* max concurrent cxx processes in the generator compile, bounded by both cores and free memory
   (peak RAM ~ jobs x per-unit RSS; a cores-only cap can OOM on a dense flow). NT_GEN_JOBS pins the
   count. Delayed (:=) so memory is read at compile time, not at package load. *)

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

(* C++ compiler for the build-time generator (g++ and clang++ take the same flags). Order: NT_GEN_CXX
   verbatim, else clang++ when present (~2.5x faster on the -O0 net-builder units), else the
   CCompilerDriver default, else "g++" on PATH. clang objects link fine against the g++-built
   libNumTracer.a (shared libstdc++ ABI on Linux). *)

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

(* tcmalloc preload for the generator RUN: the trace engine does heavy tiny-allocation churn, and
   tcmalloc_minimal is faster there and fragments less than glibc's arenas. Empty when the library
   is absent, which is also the opt-out. *)
ntTcmallocPrefix[] :=
  With[{lib = SelectFirst[
      {"/usr/lib/libtcmalloc_minimal.so", "/usr/lib64/libtcmalloc_minimal.so",
       "/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4", "/usr/lib/libtcmalloc_minimal.so.4"},
      FileExistsQ, None]},
    If[lib === None, "", "LD_PRELOAD='" <> lib <> "' "]];

(* The single definition of "does this emission target device code"; it enables gen.hpp's size-gated
   `__noinline__` (device-only). Passed explicitly, not sniffed from the generator's decorator,
   because ntKokkosDecor has already rewritten `__device__` to the Kokkos macros. `Automatic` checks
   the raw CUDA spelling; DiFfRG_compat passes True/False. Used for both the online env prefix and
   the manifest's "device" field (offline), so the two paths cannot disagree. *)
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

(* Resolve symlinks so two spellings of the same directory compare equal (else a symlinked home
   gives paths with no common prefix and the manifest's relative paths break). AbsoluteFileName needs
   the path to EXIST, and the traces header may not exist yet, so canonicalise the deepest existing
   ANCESTOR and re-append the segments below it. *)
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
   Written beside the kernel headers and COMMITTED with them. Two jobs: the "generated" 0|1 switch
   (emitting sets 0, the numtrace target sets 1 once built; a fresh clone has kernels but no gen/,
   so timestamps cannot say this) and the metadata the offline build needs. One writer per file, so
   parallel numtrace jobs never race. Paths are basenames relative to "gen_dir", itself relative to
   the flow dir. *)

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
    (* thread caps from SetNumTracerThreads at emit time (0 = unset). The offline cmake -P step
       inherits no environment, so they must be recorded here; the driver re-applies them, only
       ever lowering the build's -jN. *)
    "maxw"         -> ntEnvPosInt["NT_GEN_MAXW"],
    "maxw_b"        -> ntEnvPosInt["NT_GEN_MAXW_B"],
    (* device target: offline, this field is the only way the generator learns it (online it comes
       from ntDeviceEnvPrefix); without it the device-only noinline gate never applies *)
    "device"       -> ntDeviceTargetQ[deviceTarget, decor],
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

(* The generator compiles gen_*.cpp against the CURRENT headers but LINKS the prebuilt archive. A
   header edit without a library rebuild is a silent ODR/ABI mismatch: the generator runs, kernels.hh
   is plausible but wrong, and kernel.hh (which only names the traces) stays byte-identical. mtime is
   a coarse proxy, but the failure is silent and total, so an archive older than any header aborts. *)

MakeNTKernel::stalelib = "The engine archive\n  `1`\nis older than the header\n  `2`\nthe generator compiles against. Linking them mixes two versions of the engine types and silently produces wrong traces. Rebuild it first:\n  cmake --build <repo>/numtracer/build --target NumTracer\n(NT_GEN_LIB selects a specific archive; NT_ALLOW_STALE_LIB=1 overrides this check.)";

genLibStaleQ[lib_, incDir_] := Module[{hdrs, newest},
    hdrs = FileNames["*.hpp" | "*.h", incDir, Infinity];
    If[hdrs === {},
      Return[False]];
    newest = Max[AbsoluteTime /@ (FileDate[#, "Modification"]& /@ hdrs)];
    AbsoluteTime[FileDate[lib, "Modification"]] < newest];

(* Locate the compiled engine library libNumTracer.a for the generator LINK, so the generator only
   parses declaration headers instead of compiling the whole engine every run. Searched: NT_GEN_LIB
   (a full path), the installed <prefix>/lib, the in-tree <repo>/numtracer/build. $Failed makes the
   generator fall back to a slower header-only compile (-DNUMTRACER_HEADER_ONLY=1). *)
resolveGenLib[incDir_] := Module[{env, base, cands, lib},
    env = Environment["NT_GEN_LIB"];
    lib =
      If[StringQ[env] && FileExistsQ[env],
        env,
        (* parent of include/: <prefix> (installed) or <repo>/numtracer (in-tree) *)
        base = DirectoryName[incDir];
        cands =
          {
            (* installed layout *)
            FileNameJoin[{base, "lib", "libNumTracer.a"}],
            (* in-tree default build dir *)
            FileNameJoin[{base, "build", "libNumTracer.a"}]
          };
        SelectFirst[cands, FileExistsQ, $Failed]];
    If[lib =!= $Failed && !ntEnvFlag["NT_ALLOW_STALE_LIB"] && genLibStaleQ[lib, incDir],
      Module[{hdrs = FileNames["*.hpp" | "*.h", incDir, Infinity], newest},
        newest = First @ SortBy[hdrs, -AbsoluteTime[FileDate[#, "Modification"]]&];
        Message[MakeNTKernel::stalelib, lib, newest];
        Abort[]]];
    lib];

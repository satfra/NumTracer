(* Print the Wolfram user Applications directory, queried by CMake at configure
   time to decide where the NumTracer Mathematica package installs (same
   mechanism as FunKit's utils/get_wolfram_app_dir.m). *)
FileNameJoin[{$UserBaseDirectory, "Applications"}] // Print;
Exit[]

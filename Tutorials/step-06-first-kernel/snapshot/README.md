The headers `first_kernel.wls` generates, as quoted by the step-06 page. The tutorial tests
(`first_kernel_snapshot_*`) check that a fresh generation still produces exactly these files; after
an intended codegen change, refresh them with

    cp <build>/gen/first_kernel/first_kernel*.hh Tutorials/step-06-first-kernel/snapshot/

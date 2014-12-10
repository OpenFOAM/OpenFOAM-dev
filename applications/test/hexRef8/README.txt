Test application for volField and surfaceField mapping with
refinement/unrefinement

Run

    block/Allrun

to compile and map a few fields. Note that hexRef8 cannot be
run in inflation mode - there is the problem of getting
a set of faces to sweep that is consistent for a cell and
all its neighbours.

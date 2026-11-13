# AGENTS.md

Guidance for coding agents working in this repository.

## What this is

cfMesh for the OpenFOAM of the OpenFOAM Foundation (OpenFOAM-14), maintained
by LEC GmbH. It starts from https://gitlab.com/openfoam/community/cfmesh;
`README` lists the differences.

## Layout

- `cfMeshLibrary/` - the library: `meshGenerators/` (cartesian, cartesian2D,
  tet, voronoi) and `utilities/` (octree, surface tools, boundary layers,
  smoothers).
- `utilities/cfMesh/` - the one executable that runs all generators;
  the other folders of `utilities/` are surface and mesh tools.
- `bin/` - the scripts `cartesianMesh`, `cartesian2DMesh`, `tetMesh` and
  `pMesh`, which call `cfMesh -mode`.
- `tutorials/` - one folder per generator, each case with `Allrun`.
- `userGuide/` - the user guide.

## Build and test

With the OpenFOAM-14 environment set:

    ./Allwmake

The library and the utilities install to `$FOAM_SITE_LIBBIN` and
`$FOAM_SITE_APPBIN`. `wmake` in `cfMeshLibrary` rebuilds the library alone.

There are no unit tests. Run the tutorial cases that exercise the changed
code (`tutorials/Allrun` runs all of them) and check the result with
`checkMesh`.

## Conventions

- OpenFOAM coding style: 4 spaces, no tabs, lines of at most 80 characters.
- Parallelisation is OpenMP only; there is no MPI layer and none is added.
- Shared data written inside a parallel region needs a reason why it is
  safe: one writer per entry, `critical`, `atomic` or a reduction.
- OpenFOAM's containers and tools are used where they exist; no new
  containers and no utilities that duplicate an OpenFOAM tool.
- `meshDict` follows `snappyHexMeshDict` where the two overlap: the stage
  switches `castellatedMesh`, `snap` and `addLayers`, dictionaries rather
  than keyword and dictionary forms of one setting.
- Comments state a reason that the code does not show; they do not describe
  the history of a change.
- One change per commit; the body says what was wrong and what the code
  does now.
- A change in behaviour or a new keyword of `meshDict` is described in
  `README`.

# cfMesh for Foundation OpenFOAM

cfMesh is a library for mesh generation. This repository is cfMesh for
Foundation OpenFOAM (the OpenFOAM of the OpenFOAM Foundation, version 14),
maintained by LEC GmbH. It starts from
[cfMesh](https://gitlab.com/openfoam/community/cfmesh).

[`README.original`](README.original) is the README of the upstream cfMesh,
unchanged. Its download and contact addresses are those of upstream, not of
this repository.

## Differences to upstream

| Area | Difference |
|---|---|
| Parallelisation | Shared-memory (OpenMP) only; the MPI layer is removed. |
| Executable | One executable, `cfMesh`, runs all generators. The generator is selected with `-mode` or `meshMode` in `meshDict` (`cartesian`, `cartesian2D`, `tet`, `polyhedral`); the scripts in `bin` (`cartesianMesh`, `cartesian2DMesh`, `tetMesh`, `pMesh`) call it. |
| Stages | The stages follow `snappyHexMesh`: `castellatedMesh`, `snap` and `addLayers` are switches in `meshDict`, and the mesh is written after each stage. `-region`, `-mesh`, `-dict`, `-overwrite` and `-noOverwrite` are supported. |
| `insidePoints` | Selects which disconnected regions of the surface are meshed ([`twoBoxes`](tutorials/cartesianMesh/twoBoxes), [`threeBoxes`](tutorials/cartesianMesh/threeBoxes)). |
| `objectRefinements`, `type surface` | Refines the cells inside a closed surface. |
| `localRefinement`, `distLevels` and `distCellSizes` | Several refinement shells around one patch. |
| `meshDict` | Settings that had a keyword and a dictionary form are dictionaries only. |
| Code | cfMesh's own containers are replaced by those of Foundation OpenFOAM, and utilities that duplicate its tools are removed. |

## Building

With the environment of Foundation OpenFOAM-14 set:

```sh
./Allwmake
```

The library and the utilities install to `$FOAM_SITE_LIBBIN` and
`$FOAM_SITE_APPBIN`. Add `bin` to `PATH` for the scripts.

## Getting started

Every tutorial case has an `Allrun` and an `Allclean` script. The first
example is a simplified car body:

```sh
cd tutorials/cartesianMesh/asmoOctree
./Allrun
```

The user guide of the original cfMesh is in [`userGuide`](userGuide).

## Bug reports and contributions

- **Bug reports** go to the [issue tracker](../../issues) of this repository.
  A report that can be reproduced has the Foundation OpenFOAM version, `meshDict`, the
  surface file and the log of the mesher.
- **Pull requests** are welcome. [`AGENTS.md`](AGENTS.md) describes the
  layout of the code and the conventions for code and commits; they hold for
  every contributor.

## Licence and credits

cfMesh is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version. See [`COPYING`](COPYING).

cfMesh is the work of Franjo Juretic, Philippose Rajan and Ivor Clifford
(Creative Fields, Ltd.).

The work of LEC GmbH in this repository was done with the support of an AI
coding assistant (Claude Code by Anthropic). `AGENTS.md` holds the guidance
for such assistants.

This offering is not approved or endorsed by the OpenFOAM Foundation or by
OpenCFD Limited, the producer of the OpenFOAM software and owner of the
OPENFOAM® trade mark.

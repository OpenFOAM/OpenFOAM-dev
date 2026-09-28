/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "patchRegionSplit.H"
#include "PatchTools.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchRegionSplit, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchRegionSplit::patchRegionSplit
(
    const polyMesh& mesh,
    const labelList& faces
)
:
    regionSplitBase(faces.size())
{
    const uindirectPrimitivePatch patch
    (
        UIndirectList<face>(mesh.faces(), faces),
        mesh.points()
    );

    const label nLocalZones = PatchTools::markZones(patch, boolList(), *this);

    nRegions_ =
        Pstream::parRun()
      ? compactGlobalRegionSplit(globalIndex(nLocalZones), *this)
      : compactLocalRegionSplit(*this);
}


Foam::patchRegionSplit::patchRegionSplit(const polyPatch& patch)
:
    regionSplitBase(patch.size())
{
    const label nLocalZones = PatchTools::markZones(patch, boolList(), *this);

    nRegions_ =
        Pstream::parRun()
      ? compactGlobalRegionSplit(globalIndex(nLocalZones), *this)
      : compactLocalRegionSplit(*this);
}


// ************************************************************************* //

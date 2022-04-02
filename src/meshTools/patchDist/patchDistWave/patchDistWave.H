/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Class
    Foam::patchDistWave

Description
    Takes a set of patches to start FaceCellWave from and computed the distance
    at patches and possibly additional transported data.

SourceFiles
    patchDistWave.C

\*---------------------------------------------------------------------------*/

#ifndef patchDistWave_H
#define patchDistWave_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistWave
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Get initial set of changed faces
labelList getChangedFaces
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs
);

//- Wave distance data from faces
label wave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    scalarField& cellDistance
);

//- Calculate distance data from patches
label calculate
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistWave
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

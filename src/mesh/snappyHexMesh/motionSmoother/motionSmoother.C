/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "motionSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSmoother, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSmoother::motionSmoother
(
    polyMesh& mesh,
    pointMesh& pMesh,
    indirectPrimitivePatch& pp,
    const labelList& adaptPatchIDs,
    const dictionary& paramDict
)
:
    motionSmootherData(pMesh),
    motionSmootherAlgo
    (
        mesh,
        pMesh,
        pp,
        motionSmootherData::displacement_,
        motionSmootherData::scale_,
        motionSmootherData::oldPoints_,
        adaptPatchIDs,
        paramDict
    )
{}


Foam::motionSmoother::motionSmoother
(
    polyMesh& mesh,
    indirectPrimitivePatch& pp,
    const labelList& adaptPatchIDs,
    const pointVectorField& displacement,
    const dictionary& paramDict
)
:
    motionSmootherData(displacement),
    motionSmootherAlgo
    (
        mesh,
        const_cast<pointMesh&>(displacement.mesh()),
        pp,
        motionSmootherData::displacement_,
        motionSmootherData::scale_,
        motionSmootherData::oldPoints_,
        adaptPatchIDs,
        paramDict
    )
{}


// ************************************************************************* //

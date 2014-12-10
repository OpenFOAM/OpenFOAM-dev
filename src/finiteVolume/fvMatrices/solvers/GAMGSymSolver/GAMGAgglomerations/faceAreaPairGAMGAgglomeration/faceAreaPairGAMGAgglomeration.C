/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "faceAreaPairGAMGAgglomeration.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceAreaPairGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        faceAreaPairGAMGAgglomeration,
        lduMesh
    );

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        faceAreaPairGAMGAgglomeration,
        geometry
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceAreaPairGAMGAgglomeration::faceAreaPairGAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(mesh, controlDict)
{
    const fvMesh& fvmesh = refCast<const fvMesh>(mesh);

    //agglomerate(mesh, sqrt(fvmesh.magSf().internalField()));
    agglomerate
    (
        mesh,
        mag
        (
            cmptMultiply
            (
                fvmesh.Sf().internalField()
               /sqrt(fvmesh.magSf().internalField()),
                vector(1, 1.01, 1.02)
                //vector::one
            )
        )
    );
}


Foam::faceAreaPairGAMGAgglomeration::faceAreaPairGAMGAgglomeration
(
    const lduMesh& mesh,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(mesh, controlDict)
{
    //agglomerate(mesh, sqrt(mag(faceAreas)));
    agglomerate
    (
        mesh,
        mag
        (
            cmptMultiply
            (
                faceAreas
               /sqrt(mag(faceAreas)),
                vector(1, 1.01, 1.02)
                //vector::one
            )
        )
    );
}


// ************************************************************************* //

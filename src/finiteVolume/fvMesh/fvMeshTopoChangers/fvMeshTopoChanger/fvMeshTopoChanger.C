/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "fvMeshTopoChanger.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshTopoChanger, 0);
    defineRunTimeSelectionTable(fvMeshTopoChanger, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChanger::fvMeshTopoChanger(fvMesh& mesh)
:
    mesh_(mesh),
    dynamicMeshDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh.dbDir(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
    )
{}


Foam::fvMeshTopoChanger::velocityMotionCorrection::velocityMotionCorrection
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    velocityFields_(dict.lookupOrDefault("velocityFields", wordList()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChanger::~fvMeshTopoChanger()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshTopoChanger::velocityMotionCorrection::update() const
{
    forAll(velocityFields_, i)
    {
        if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
        {
            mesh_.lookupObjectRef<volVectorField>
            (
                velocityFields_[i]
            ).correctBoundaryConditions();
        }
    }
}

// ************************************************************************* //

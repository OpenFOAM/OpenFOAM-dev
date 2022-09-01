/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "buoyancy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::buoyancy::buoyancy(const fvMesh& mesh_)
:
    mesh(mesh_),

    runTime(mesh.time()),

    g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, 0)
    ),

    pRef
    (
        IOobject
        (
            "pRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimPressure, 0)
    ),

    ghRef(-mag(g)*hRef),

    gh("gh", (g & mesh.C()) - ghRef),

    ghf("ghf", (g & mesh.Cf()) - ghRef),

    p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    mesh.schemes().setFluxRequired(p_rgh.name());
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solvers::buoyancy> Foam::solvers::buoyancy::New
(
    const fvMesh& mesh
)
{
    return typeIOobject<volScalarField>
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        ).headerOk()
      ? autoPtr<buoyancy>(new solvers::buoyancy(mesh))
      : autoPtr<buoyancy>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::buoyancy::~buoyancy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::buoyancy::moveMesh()
{
    gh = (g & mesh.C()) - ghRef;
    ghf = (g & mesh.Cf()) - ghRef;
}


// ************************************************************************* //

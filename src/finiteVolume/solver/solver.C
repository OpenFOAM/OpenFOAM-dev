/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "solver.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solver, 0);
    defineRunTimeSelectionTable(solver, fvMesh);
}


Foam::scalar Foam::solver::deltaTFactor = 1.2;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solver::writeData(Ostream&) const
{
    NotImplemented;
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solver::solver(fvMesh& mesh_)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh_.time().timeName(),
            mesh_
        )
    ),

    mesh(mesh_),

    runTime(mesh.time()),

    pimple(mesh),

    steady(mesh.schemes().steady()),

    LTS(fv::localEulerDdt::enabled(mesh)),

    fvModelsPtr(nullptr),
    fvConstraintsPtr(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solver::~solver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::fvModels& Foam::solver::fvModels() const
{
    if (!fvModelsPtr)
    {
        fvModelsPtr = &Foam::fvModels::New(mesh);
    }

    return *fvModelsPtr;
}


Foam::fvConstraints& Foam::solver::fvConstraints() const
{
    if (!fvConstraintsPtr)
    {
        fvConstraintsPtr = &Foam::fvConstraints::New(mesh);
    }

    return *fvConstraintsPtr;
}


// ************************************************************************* //

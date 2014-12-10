/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "wallDependentModel.H"
#include "wallDist.H"
#include "wallDistReflection.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDependentModel::wallDependentModel(const fvMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDependentModel::~wallDependentModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::wallDependentModel::yWall() const
{
    if (!mesh_.foundObject<volScalarField>("yWall"))
    {
        wallDist w(mesh_);

        volScalarField* yPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "yWall",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                w.y()
            )
        );

        yPtr->checkIn();
    }

    return mesh_.lookupObject<volScalarField>("yWall");
}


const Foam::volVectorField& Foam::wallDependentModel::nWall() const
{
    if (!mesh_.foundObject<volVectorField>("nWall"))
    {
        wallDistReflection w(mesh_);

        if (!mesh_.foundObject<volScalarField>("yWall"))
        {
            volScalarField* yPtr
            (
                new volScalarField
                (
                    IOobject
                    (
                        "yWall",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        true
                    ),
                    w.y()
                )
            );

            yPtr->checkIn();
        }

        volVectorField* nPtr
        (
            new volVectorField
            (
                IOobject
                (
                    "nWall",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                w.n()
            )
        );

        nPtr->checkIn();
    }

    return mesh_.lookupObject<volVectorField>("nWall");
}


// ************************************************************************* //

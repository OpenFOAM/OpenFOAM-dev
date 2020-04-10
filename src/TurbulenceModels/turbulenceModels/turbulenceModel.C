/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "turbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(turbulenceModel, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::IOdictionary Foam::turbulenceModel::readModelDict
(
    const objectRegistry& obr,
    const word& group,
    bool registerObject
)
{
    IOobject momentumTransport
    (
        IOobject::groupName(typeName, group),
        obr.time().constant(),
        obr,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        registerObject
    );

    if (momentumTransport.typeHeaderOk<IOdictionary>(true))
    {
        return momentumTransport;
    }
    else
    {
        IOobject momentumTransport
        (
            IOobject::groupName("momentumTransport", group),
            obr.time().constant(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            registerObject
        );

        if (momentumTransport.typeHeaderOk<IOdictionary>(true))
        {
            return momentumTransport;
        }
        else
        {
            return momentumTransport;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulenceModel::turbulenceModel
(
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi
)
:
    IOdictionary(readModelDict(U.db(), alphaRhoPhi.group(), true)),

    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    alphaRhoPhi_(alphaRhoPhi),
    phi_(phi),
    y_(mesh_)
{
    // Ensure name of IOdictionary is typeName
    rename(IOobject::groupName(typeName, alphaRhoPhi.group()));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::turbulenceModel::phi() const
{
    return phi_;
}


bool Foam::turbulenceModel::read()
{
    return regIOobject::read();
}


void Foam::turbulenceModel::validate()
{}


void Foam::turbulenceModel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


// ************************************************************************* //

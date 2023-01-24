/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "momentumTransportModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallFvPatch.H"
#include "nearWallDist.H"
#include "fvcFlux.H"
#include "fvmDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentumTransportModel, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::typeIOobject<Foam::IOdictionary>
Foam::momentumTransportModel::readModelDict
(
    const objectRegistry& obr,
    const word& group,
    bool registerObject
)
{
    typeIOobject<IOdictionary> momentumTransport
    (
        IOobject::groupName(typeName, group),
        obr.time().constant(),
        obr,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        registerObject
    );

    if (momentumTransport.headerOk())
    {
        return momentumTransport;
    }
    else
    {
        typeIOobject<IOdictionary> turbulenceProperties
        (
            IOobject::groupName("turbulenceProperties", group),
            obr.time().constant(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            registerObject
        );

        if (turbulenceProperties.headerOk())
        {
            return turbulenceProperties;
        }
        else
        {
            return momentumTransport;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentumTransportModel::momentumTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    IOdictionary(readModelDict(U.db(), alphaRhoPhi.group(), true)),

    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    alphaRhoPhi_(alphaRhoPhi),
    phi_(phi),
    viscosity_(viscosity)
{
    // Ensure name of IOdictionary is typeName
    rename(IOobject::groupName(typeName, alphaRhoPhi.group()));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::momentumTransportModel::phi() const
{
    return phi_;
}


const Foam::volScalarField::Boundary& Foam::momentumTransportModel::y() const
{
    return nearWallDist::New(mesh_).y();
}


bool Foam::momentumTransportModel::read()
{
    return regIOobject::read();
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::momentumTransportModel::divDevTauCorr
(
    const tmp<volTensorField>& devTauCorr,
    volVectorField& U
) const
{
    return fvm::divc
    (
        tmp<surfaceVectorField>
        (
            new surfaceVectorField
            (
                groupName("devTauCorrFlux"),
                fvc::flux(devTauCorr)
            )
        ),
        U
    );
}


void Foam::momentumTransportModel::validate()
{}


void Foam::momentumTransportModel::predict()
{}


void Foam::momentumTransportModel::correct()
{}


// ************************************************************************* //

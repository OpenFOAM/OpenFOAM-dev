/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "specieFlux.H"
#include "fluidThermophysicalTransportModel.H"
#include "fvcFlux.H"
#include "multicomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(specieFlux, 0);
    addToRunTimeSelectionTable(functionObject, specieFlux, dictionary);

    defineTypeNameAndDebug(specieAdvectiveFlux, 0);
    addToRunTimeSelectionTable(functionObject, specieAdvectiveFlux, dictionary);

    defineTypeNameAndDebug(specieDiffusiveFlux, 0);
    addToRunTimeSelectionTable(functionObject, specieDiffusiveFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::specieFluxBase::calc()
{
    const word thermoName =
        IOobject::groupName
        (
            physicalProperties::typeName,
            IOobject::group(fieldName_)
        );

    const word ttmName =
        IOobject::groupName
        (
            fluidThermophysicalTransportModel::typeName,
            IOobject::group(fieldName_)
        );

    if
    (
        !foundObject<multicomponentThermo>(thermoName)
     || !foundObject<fluidThermophysicalTransportModel>(ttmName)
    )
    {
        return false;
    }

    const multicomponentThermo& thermo =
        lookupObject<multicomponentThermo>(thermoName);

    const fluidThermophysicalTransportModel& ttm =
        lookupObject<fluidThermophysicalTransportModel>(ttmName);

    const volScalarField& Yi = thermo.Y(fieldName_);

    return store(resultName_, calc(ttm, Yi));
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieFlux::calc
(
    const fluidThermophysicalTransportModel& ttm,
    const volScalarField& Yi
)
{
    return calcPhiYif(ttm, Yi) + calcJ(ttm, Yi);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieAdvectiveFlux::calc
(
    const fluidThermophysicalTransportModel& ttm,
    const volScalarField& Yi
)
{
    return calcPhiYif(ttm, Yi);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieDiffusiveFlux::calc
(
    const fluidThermophysicalTransportModel& ttm,
    const volScalarField& Yi
)
{
    return calcJ(ttm, Yi);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieFluxBase::calcPhiYif
(
    const fluidThermophysicalTransportModel& ttm,
    const volScalarField& Yi
) const
{
    const surfaceScalarField& phi = ttm.momentumTransport().alphaRhoPhi();

    return
        fvc::flux
        (
            phi,
            Yi,
            "div(" + phi.name() + "," + schemesField_ + ")"
        )/Yi.mesh().magSf();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieFluxBase::calcJ
(
    const fluidThermophysicalTransportModel& ttm,

    const volScalarField& Yi
) const
{
    return ttm.j(Yi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::specieFluxBase::specieFluxBase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& typeName
)
:
    fieldExpression(name, runTime, dict, typeName),
    schemesField_(dict.lookupOrDefault<word>("schemesField", "Yi"))
{}


Foam::functionObjects::specieFlux::specieFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    specieFluxBase(name, runTime, dict, typeName)
{}


Foam::functionObjects::specieAdvectiveFlux::specieAdvectiveFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    specieFluxBase(name, runTime, dict, typeName)
{}


Foam::functionObjects::specieDiffusiveFlux::specieDiffusiveFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    specieFluxBase(name, runTime, dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::specieFluxBase::~specieFluxBase()
{}


// ************************************************************************* //

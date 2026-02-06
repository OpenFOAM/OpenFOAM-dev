/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "energyFlux.H"
#include "fluidThermophysicalTransportModel.H"
#include "fvcFlux.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(energyFlux, 0);
    addToRunTimeSelectionTable(functionObject, energyFlux, dictionary);

    defineTypeNameAndDebug(energyAdvectiveFlux, 0);
    addToRunTimeSelectionTable(functionObject, energyAdvectiveFlux, dictionary);

    defineTypeNameAndDebug(heatFlux, 0);
    addToRunTimeSelectionTable(functionObject, heatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::energyFluxBase::calc()
{
    const word ttmName =
        IOobject::groupName
        (
            thermophysicalTransportModel::typeName,
            phaseName_
        );

    if (!foundObject<thermophysicalTransportModel>(ttmName))
    {
        return false;
    }

    const thermophysicalTransportModel& ttm =
        lookupObject<thermophysicalTransportModel>(ttmName);

    store(resultName_, calc(ttm));

    return true;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::energyFlux::calc
(
    const thermophysicalTransportModel& ttm
) const
{
    return calcPhihs(ttm) + calcQ(ttm);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::energyAdvectiveFlux::calc
(
    const thermophysicalTransportModel& ttm
) const
{
    return calcPhihs(ttm);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::heatFlux::calc
(
    const thermophysicalTransportModel& ttm
) const
{
    return calcQ(ttm);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::energyFluxBase::calcPhihs
(
    const thermophysicalTransportModel& ttm
) const
{
    if (isA<fluidThermophysicalTransportModel>(ttm))
    {
        const fluidThermophysicalTransportModel& fttm =
            refCast<const fluidThermophysicalTransportModel>(ttm);

        const surfaceScalarField& phi = fttm.momentumTransport().alphaRhoPhi();

        const word& schemesField =
            schemesField_ != word::null
          ? schemesField_
          : fttm.thermo().he().name();

        return
            fvc::flux
            (
                phi,
                fttm.thermo().hs(),
                "div(" + phi.name() + "," + schemesField + ")"
            );
    }

    return surfaceScalarField::New("0", mesh(), dimPower);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::energyFluxBase::calcQ
(
    const thermophysicalTransportModel& ttm
) const
{
    return mesh().magSf()*ttm.q();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::energyFluxBase::energyFluxBase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& typeName
)
:
    fieldExpression(name, runTime, dict, typeName, noFieldName_),
    schemesField_(dict.lookupOrDefault<word>("schemesField", word::null)),
    phaseName_(word::null)
{
    read(dict);
}


Foam::functionObjects::energyFlux::energyFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    energyFluxBase(name, runTime, dict, typeName)
{}


Foam::functionObjects::energyAdvectiveFlux::energyAdvectiveFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    energyFluxBase(name, runTime, dict, typeName)
{}


Foam::functionObjects::heatFlux::heatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    energyFluxBase(name, runTime, dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::energyFluxBase::~energyFluxBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::energyFluxBase::read
(
    const dictionary& dict
)
{
    fieldExpression::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    return true;
}


// ************************************************************************* //

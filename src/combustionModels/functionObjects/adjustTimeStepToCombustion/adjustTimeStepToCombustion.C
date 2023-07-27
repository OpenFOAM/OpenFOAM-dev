/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "adjustTimeStepToCombustion.H"
#include "combustionModel.H"
#include "solver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(adjustTimeStepToCombustion, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        adjustTimeStepToCombustion,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::typeIOobject<Foam::timeIOdictionary>
Foam::functionObjects::adjustTimeStepToCombustion::propsDictIo
(
    const IOobject::readOption& r
) const
{
    return
        typeIOobject<timeIOdictionary>
        (
            name() + "Properties",
            obr_.time().name(),
            "uniform",
            obr_,
            r,
            IOobject::NO_WRITE,
            false
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::adjustTimeStepToCombustion::adjustTimeStepToCombustion
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    phaseName_(word::null),
    maxCo_(NaN),
    extrapolate_(false),
    haveCombustionDeltaT0_(false),
    combustionDeltaT0_(NaN)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::adjustTimeStepToCombustion::~adjustTimeStepToCombustion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::adjustTimeStepToCombustion::read
(
    const dictionary& dict
)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);
    maxCo_ = dict.lookupOrDefault<scalar>("maxCo", 1);
    extrapolate_ = dict.lookupOrDefault<bool>("extrapolate", false);

    typeIOobject<timeIOdictionary> propsDictIo
    (
        this->propsDictIo(IOobject::MUST_READ_IF_MODIFIED)
    );

    if (propsDictIo.headerOk())
    {
        const timeIOdictionary propsDict(propsDictIo);

        haveCombustionDeltaT0_ = true;
        combustionDeltaT0_ = propsDict.lookup<scalar>("combustionDeltaT");
    }
    else
    {
        haveCombustionDeltaT0_ = false;
    }

    return true;
}


bool Foam::functionObjects::adjustTimeStepToCombustion::execute()
{
    return true;
}


bool Foam::functionObjects::adjustTimeStepToCombustion::write()
{
    if (extrapolate_ && obr_.time().writeTime())
    {
        timeIOdictionary propsDict(propsDictIo(IOobject::NO_READ));

        propsDict.add("combustionDeltaT", combustionDeltaT0_);

        propsDict.regIOobject::write();
    }

    return true;
}


Foam::scalar
Foam::functionObjects::adjustTimeStepToCombustion::maxDeltaT() const
{
    if (!time_.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        return vGreat;
    }

    const combustionModel& combustion =
        obr_.lookupObject<combustionModel>
        (
            IOobject::groupName
            (
                combustionModel::combustionPropertiesName,
                phaseName_
            )
        );

    const fluidMulticomponentThermo& thermo = combustion.thermo();

    // Build a mass turnover rate
    volScalarField::Internal rhoDotByRho
    (
        volScalarField::Internal::New
        (
            "rhoDotByRho",
            combustion.mesh(),
            dimensionedScalar(dimless/dimTime, 0)
        )
    );
    forAll(thermo.Y(), i)
    {
        if (thermo.solveSpecie(i))
        {
            rhoDotByRho += mag(combustion.R(i))/2/thermo.rho()();
        }
    }

    // Convert to a time-scale
    const scalar combustionDeltaT1 = maxCo_/max(gMax(rhoDotByRho), vSmall);

    // We want to clip the time-step to the time-scale, but also additionally
    // reduce the time-step significantly if that time-scale is reducing
    // rapidly. This helps us catch the onset of reactions. The following does
    // this by passing the change ratio into a steeply non-linear function.
    scalar deltaT;
    if (extrapolate_ && haveCombustionDeltaT0_)
    {
        const scalar x = min(max(combustionDeltaT1/combustionDeltaT0_, 0), 1);
        const scalar f = (1 - rootSmall)*(1 - sqrt(1 - sqr(x))) + rootSmall;
        deltaT = f*combustionDeltaT1;
    }
    else
    {
        deltaT = combustionDeltaT1;
    }

    // Store the latest combustion time-step
    haveCombustionDeltaT0_ = true;
    combustionDeltaT0_ = combustionDeltaT1;

    return deltaT;
}


// ************************************************************************* //

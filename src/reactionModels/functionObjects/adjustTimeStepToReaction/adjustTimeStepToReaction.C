/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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

#include "adjustTimeStepToReaction.H"
#include "reactionModel.H"
#include "solver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(adjustTimeStepToReaction, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        adjustTimeStepToReaction,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::typeIOobject<Foam::timeIOdictionary>
Foam::functionObjects::adjustTimeStepToReaction::propsDictIo
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

Foam::functionObjects::adjustTimeStepToReaction::adjustTimeStepToReaction
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
    haveReactionDeltaT0_(false),
    reactionDeltaT0_(NaN)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::adjustTimeStepToReaction::~adjustTimeStepToReaction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::adjustTimeStepToReaction::read
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

        haveReactionDeltaT0_ = true;
        reactionDeltaT0_ = propsDict.lookup<scalar>("reactionDeltaT");
    }
    else
    {
        haveReactionDeltaT0_ = false;
    }

    return true;
}


bool Foam::functionObjects::adjustTimeStepToReaction::execute()
{
    return true;
}


bool Foam::functionObjects::adjustTimeStepToReaction::write()
{
    if (extrapolate_ && obr_.time().writeTime())
    {
        timeIOdictionary propsDict(propsDictIo(IOobject::NO_READ));

        propsDict.add("reactionDeltaT", reactionDeltaT0_);

        propsDict.regIOobject::write();
    }

    return true;
}


Foam::scalar
Foam::functionObjects::adjustTimeStepToReaction::maxDeltaT() const
{
    if (!time_.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        return vGreat;
    }

    const reactionModel& reaction =
        obr_.lookupObject<reactionModel>
        (
            IOobject::groupName
            (
                reactionModel::reactionPropertiesName,
                phaseName_
            )
        );

    const fluidMulticomponentThermo& thermo = reaction.thermo();

    // Build a mass turnover rate
    volScalarField::Internal rhoDotByRho
    (
        volScalarField::Internal::New
        (
            "rhoDotByRho",
            reaction.mesh(),
            dimensionedScalar(dimless/dimTime, 0)
        )
    );
    forAll(thermo.Y(), i)
    {
        if (thermo.solveSpecie(i))
        {
            rhoDotByRho += mag(reaction.R(i))/2/thermo.rho()();
        }
    }

    // Convert to a time-scale
    const scalar reactionDeltaT1 = maxCo_/max(gMax(rhoDotByRho), vSmall);

    // We want to clip the time-step to the time-scale, but also additionally
    // reduce the time-step significantly if that time-scale is reducing
    // rapidly. This helps us catch the onset of reactions. The following does
    // this by passing the change ratio into a steeply non-linear function.
    scalar deltaT;
    if (extrapolate_ && haveReactionDeltaT0_)
    {
        const scalar x = min(max(reactionDeltaT1/reactionDeltaT0_, 0), 1);
        const scalar f = (1 - rootSmall)*(1 - sqrt(1 - sqr(x))) + rootSmall;
        deltaT = f*reactionDeltaT1;
    }
    else
    {
        deltaT = reactionDeltaT1;
    }

    // Store the latest reaction time-step
    haveReactionDeltaT0_ = true;
    reactionDeltaT0_ = reactionDeltaT1;

    return deltaT;
}


// ************************************************************************* //

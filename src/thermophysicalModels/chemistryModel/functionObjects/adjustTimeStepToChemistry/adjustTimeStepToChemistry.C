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

#include "adjustTimeStepToChemistry.H"
#include "basicChemistryModel.H"
#include "solver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(adjustTimeStepToChemistry, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        adjustTimeStepToChemistry,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::adjustTimeStepToChemistry::adjustTimeStepToChemistry
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    phaseName_(word::null)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::adjustTimeStepToChemistry::~adjustTimeStepToChemistry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::adjustTimeStepToChemistry::read
(
    const dictionary& dict
)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    return true;
}


bool Foam::functionObjects::adjustTimeStepToChemistry::execute()
{
    if (!time_.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        return true;
    }

    const basicChemistryModel& chemistry =
        obr_.lookupObject<basicChemistryModel>
        (
            IOobject::groupName("chemistryProperties", phaseName_)
        );

    const scalar deltaT = gMin(chemistry.deltaTChem());

    // The solver has not adjusted the time-step yet. When it does, if it is
    // within the physical and specified limits it will increase it by a
    // fixed factor. So, we clip it here to the chemical time-step divided by
    // that factor. The solver will then increase it to the chemical time-step
    // if it can.
    const_cast<Time&>(time_).setDeltaTNoAdjust
    (
        min(deltaT/solver::deltaTFactor, time_.deltaTValue())
    );

    return true;
}


bool Foam::functionObjects::adjustTimeStepToChemistry::write()
{
    return true;
}


// ************************************************************************* //

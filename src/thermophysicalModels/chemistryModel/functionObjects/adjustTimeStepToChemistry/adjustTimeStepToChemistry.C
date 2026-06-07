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

#include "adjustTimeStepToChemistry.H"
#include "chemistryModel.H"
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
    return true;
}


bool Foam::functionObjects::adjustTimeStepToChemistry::write()
{
    return true;
}


Foam::scalar Foam::functionObjects::adjustTimeStepToChemistry::maxDeltaT() const
{
    if (!time_.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        return vGreat;
    }

    const chemistryModel& chemistry =
        obr_.lookupObject<chemistryModel>
        (
            IOobject::groupName("chemistryProperties", phaseName_)
        );

    return gMin(chemistry.deltaTChem());
}


// ************************************************************************* //

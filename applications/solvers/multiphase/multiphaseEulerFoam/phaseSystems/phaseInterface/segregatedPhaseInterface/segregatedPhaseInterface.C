/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "segregatedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    bool segregatedPhaseInterfaceAddedHeadSeparator =
        phaseInterface::addHeadSeparator(segregatedPhaseInterface::separator());

    bool segregatedPhaseInterfaceAddedOldSeparatorToSeparator =
        phaseInterface::addOldSeparatorToSeparator
        (
            "and",
            segregatedPhaseInterface::separator()
        );
}

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        segregatedPhaseInterface,
        separatorsToTypeName({separator()}).c_str(),
        0
    );
    addToRunTimeSelectionTable(phaseInterface, segregatedPhaseInterface, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::segregatedPhaseInterface::segregatedPhaseInterface
(
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    phaseInterface(phase1, phase2)
{}


Foam::segregatedPhaseInterface::segregatedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::segregatedPhaseInterface::~segregatedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::segregatedPhaseInterface::name() const
{
    return phase1().name() + '_' + separator() + '_' + phase2().name();
}


// ************************************************************************* //

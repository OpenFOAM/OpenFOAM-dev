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

#include "segregatedSidedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        segregatedSidedPhaseInterface,
        separatorsToTypeName
        ({
            segregatedPhaseInterface::separator(),
            sidedPhaseInterface::separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        phaseInterface,
        segregatedSidedPhaseInterface,
        word
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::segregatedSidedPhaseInterface::segregatedSidedPhaseInterface
(
    const phaseModel& phase,
    const phaseModel& otherPhase
)
:
    phaseInterface(phase, otherPhase),
    segregatedPhaseInterface(phase, otherPhase),
    sidedPhaseInterface(phase, otherPhase)
{}


Foam::segregatedSidedPhaseInterface::segregatedSidedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    segregatedPhaseInterface(fluid, name),
    sidedPhaseInterface(fluid, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::segregatedSidedPhaseInterface::~segregatedSidedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::segregatedSidedPhaseInterface::name() const
{
    return
        segregatedPhaseInterface::name()
      + '_'
      + sidedPhaseInterface::separator()
      + '_'
      + phase().name();
}


// ************************************************************************* //

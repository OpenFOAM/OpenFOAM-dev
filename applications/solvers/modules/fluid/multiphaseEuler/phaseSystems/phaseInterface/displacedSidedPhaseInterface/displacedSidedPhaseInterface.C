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

#include "displacedSidedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        displacedSidedPhaseInterface,
        separatorsToTypeName
        ({
            displacedPhaseInterface::separator(),
            sidedPhaseInterface::separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        phaseInterface,
        displacedSidedPhaseInterface,
        word
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacedSidedPhaseInterface::displacedSidedPhaseInterface
(
    const phaseModel& phase,
    const phaseModel& otherPhase,
    const phaseModel& displacing
)
:
    phaseInterface(phase, otherPhase),
    displacedPhaseInterface(phase, otherPhase, displacing),
    sidedPhaseInterface(phase, otherPhase)
{}


Foam::displacedSidedPhaseInterface::displacedSidedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    displacedPhaseInterface(fluid, name),
    sidedPhaseInterface(fluid, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacedSidedPhaseInterface::~displacedSidedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::displacedSidedPhaseInterface::name() const
{
    return
        displacedPhaseInterface::name()
      + '_'
      + sidedPhaseInterface::separator()
      + '_'
      + phase().name();
}


// ************************************************************************* //

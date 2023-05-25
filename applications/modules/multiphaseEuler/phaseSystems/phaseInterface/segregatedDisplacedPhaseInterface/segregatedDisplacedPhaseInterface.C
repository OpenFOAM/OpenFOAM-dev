/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "segregatedDisplacedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        segregatedDisplacedPhaseInterface,
        separatorsToTypeName
        ({
            segregatedPhaseInterface::separator(),
            displacedPhaseInterface::separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        phaseInterface,
        segregatedDisplacedPhaseInterface,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::segregatedDisplacedPhaseInterface::segregatedDisplacedPhaseInterface
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const phaseModel& displacing
)
:
    phaseInterface(phase1, phase2),
    segregatedPhaseInterface(phase1, phase2),
    displacedPhaseInterface(phase1, phase2, displacing)
{}


Foam::segregatedDisplacedPhaseInterface::segregatedDisplacedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    segregatedPhaseInterface(fluid, name),
    displacedPhaseInterface(fluid, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::segregatedDisplacedPhaseInterface::~segregatedDisplacedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::segregatedDisplacedPhaseInterface::name() const
{
    return
        segregatedPhaseInterface::name()
      + '_'
      + displacedPhaseInterface::separator()
      + '_'
      + displacing().name();
}


// ************************************************************************* //

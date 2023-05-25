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

#include "displacedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        displacedPhaseInterface,
        separatorsToTypeName
        ({
            phaseInterface::separator(),
            separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable(phaseInterface, displacedPhaseInterface, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacedPhaseInterface::displacedPhaseInterface
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const phaseModel& displacing
)
:
    phaseInterface(phase1, phase2),
    displacing_(displacing)
{}


Foam::displacedPhaseInterface::displacedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    displacing_(identifyPhases(fluid, name, {separator()}).second())
{
    if (contains(displacing_))
    {
        FatalErrorInFunction
            << "Interface " << name << " is not valid. An interface cannot "
            << "be displaced by one of its own phases." << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacedPhaseInterface::~displacedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::displacedPhaseInterface::name() const
{
    return
        phaseInterface::name()
      + '_'
      + separator()
      + '_'
      + displacing().name();
}


const Foam::phaseModel& Foam::displacedPhaseInterface::displacing() const
{
    return displacing_;
}


// ************************************************************************* //

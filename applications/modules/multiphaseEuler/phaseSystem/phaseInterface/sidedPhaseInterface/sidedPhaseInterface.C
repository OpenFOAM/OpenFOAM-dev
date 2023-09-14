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

#include "phaseInterface.H"
#include "sidedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        sidedPhaseInterface,
        separatorsToTypeName
        ({
            phaseInterface::separator(),
            separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable(phaseInterface, sidedPhaseInterface, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sidedPhaseInterface::sidedPhaseInterface
(
    const phaseModel& phase,
    const phaseModel& otherPhase
)
:
    phaseInterface(phase, otherPhase),
    phase_(phase)
{}


Foam::sidedPhaseInterface::sidedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    phase_(identifyPhases(fluid, name, {separator()}).second())
{
    if (!contains(phase_))
    {
        FatalErrorInFunction
            << "Interface " << name << " is not valid. An interface cannot "
            << "have a side that is not one of its own phases."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sidedPhaseInterface::~sidedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::sidedPhaseInterface::name() const
{
    return phaseInterface::name() + '_' + separator() + '_' + phase().name();
}


const Foam::phaseModel& Foam::sidedPhaseInterface::phase() const
{
    return phase_;
}


const Foam::phaseModel& Foam::sidedPhaseInterface::otherPhase() const
{
    return phaseInterface::otherPhase(phase_);
}


Foam::autoPtr<Foam::phaseInterface>
Foam::sidedPhaseInterface::otherInterface() const
{
    wordList nameParts = phaseInterface::nameToNameParts(fluid(), name());

    const label i =
        findIndex(nameParts, sidedPhaseInterface::separator());

    nameParts[i+1] = otherPhase().name();

    return phaseInterface::New
    (
        fluid(),
        phaseInterface::namePartsToName(fluid(), nameParts)
    );
}


// ************************************************************************* //

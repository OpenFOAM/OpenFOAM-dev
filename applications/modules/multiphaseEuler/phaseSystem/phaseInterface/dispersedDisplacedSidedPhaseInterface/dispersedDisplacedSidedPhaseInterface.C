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

#include "dispersedDisplacedSidedPhaseInterface.H"
#include "dispersedDisplacedPhaseInterface.H"
#include "displacedSidedPhaseInterface.H"
#include "dispersedSidedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        dispersedDisplacedSidedPhaseInterface,
        separatorsToTypeName
        ({
            dispersedPhaseInterface::separator(),
            displacedPhaseInterface::separator(),
            sidedPhaseInterface::separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        phaseInterface,
        dispersedDisplacedSidedPhaseInterface,
        word
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersedDisplacedSidedPhaseInterface::
dispersedDisplacedSidedPhaseInterface
(
    const phaseModel& dispersed,
    const phaseModel& continuous,
    const phaseModel& displacing,
    const phaseModel& phase
)
:
    phaseInterface(dispersed, continuous),
    dispersedPhaseInterface(dispersed, continuous),
    displacedPhaseInterface(dispersed, continuous, displacing),
    sidedPhaseInterface(phase, phaseInterface::otherPhase(phase))
{}


Foam::dispersedDisplacedSidedPhaseInterface::
dispersedDisplacedSidedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    dispersedPhaseInterface(fluid, name),
    displacedPhaseInterface(fluid, name),
    sidedPhaseInterface(fluid, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dispersedDisplacedSidedPhaseInterface::
~dispersedDisplacedSidedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::dispersedDisplacedSidedPhaseInterface::name() const
{
    return
        dispersedPhaseInterface::name()
      + '_'
      + displacedPhaseInterface::separator()
      + '_'
      + displacing().name()
      + '_'
      + sidedPhaseInterface::separator()
      + '_'
      + phase().name();
}


// ************************************************************************* //

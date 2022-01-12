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

#include "dispersedDisplacedPhaseInterface.H"
#include "dispersedPhaseInterface.H"
#include "displacedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        dispersedDisplacedPhaseInterface,
        separatorsToTypeName
        ({
            dispersedPhaseInterface::separator(),
            displacedPhaseInterface::separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        phaseInterface,
        dispersedDisplacedPhaseInterface,
        word
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersedDisplacedPhaseInterface::dispersedDisplacedPhaseInterface
(
    const phaseModel& dispersed,
    const phaseModel& continuous,
    const phaseModel& displacing
)
:
    phaseInterface(dispersed, continuous),
    dispersedPhaseInterface(dispersed, continuous),
    displacedPhaseInterface(dispersed, continuous, displacing)
{}


Foam::dispersedDisplacedPhaseInterface::dispersedDisplacedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    dispersedPhaseInterface(fluid, name),
    displacedPhaseInterface(fluid, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dispersedDisplacedPhaseInterface::~dispersedDisplacedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::dispersedDisplacedPhaseInterface::name() const
{
    return
        dispersedPhaseInterface::name()
      + '_'
      + displacedPhaseInterface::separator()
      + '_'
      + displacing().name();
}


// ************************************************************************* //

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

#include "dispersedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    bool dispersedPhaseInterfaceAddedHeadSeparator =
        phaseInterface::addHeadSeparator(dispersedPhaseInterface::separator());

    bool dispersedPhaseInterfaceAddedOldSeparatorToSeparator =
        phaseInterface::addOldSeparatorToSeparator
        (
            "in",
            dispersedPhaseInterface::separator()
        );
}

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        dispersedPhaseInterface,
        separatorsToTypeName({separator()}).c_str(),
        0
    );
    addToRunTimeSelectionTable(phaseInterface, dispersedPhaseInterface, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersedPhaseInterface::dispersedPhaseInterface
(
    const phaseModel& dispersed,
    const phaseModel& continuous
)
:
    phaseInterface(dispersed, continuous),
    dispersed_(dispersed)
{}


Foam::dispersedPhaseInterface::dispersedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    dispersed_(identifyPhases(fluid, name, {separator()}).first())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dispersedPhaseInterface::~dispersedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::dispersedPhaseInterface::name() const
{
    return dispersed().name() + '_' + separator() + '_' + continuous().name();
}


const Foam::phaseModel& Foam::dispersedPhaseInterface::dispersed() const
{
    return dispersed_;
}


const Foam::phaseModel& Foam::dispersedPhaseInterface::continuous() const
{
    return otherPhase(dispersed_);
}


Foam::tmp<Foam::volVectorField> Foam::dispersedPhaseInterface::Ur() const
{
    return dispersed().U() - continuous().U();
}


Foam::tmp<Foam::volScalarField> Foam::dispersedPhaseInterface::Re() const
{
    return magUr()*dispersed().d()/continuous().thermo().nu();
}


Foam::tmp<Foam::volScalarField> Foam::dispersedPhaseInterface::Pr() const
{
    return
         continuous().thermo().nu()
        *continuous().thermo().Cp()
        *continuous().rho()
        /continuous().thermo().kappa();
}


Foam::tmp<Foam::volScalarField> Foam::dispersedPhaseInterface::Eo() const
{
    return Eo(dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::dispersedPhaseInterface::Eo
(
    const volScalarField& d
) const
{
    return
        mag(dispersed().rho() - continuous().rho())
       *mag(g())
       *sqr(d)
       /sigma();
}


Foam::tmp<Foam::volScalarField> Foam::dispersedPhaseInterface::Mo() const
{
    return
        mag(g())
       *continuous().thermo().nu()
       *pow3
        (
            continuous().thermo().nu()
           *continuous().rho()
           /sigma()
        );
}


Foam::tmp<Foam::volScalarField> Foam::dispersedPhaseInterface::Ta() const
{
    return Re()*pow(Mo(), 0.23);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "dispersedDragModel.H"
#include "phaseSystem.H"
#include "noSwarm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::dispersedDragModel::dispersedDragModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dragModel(dict, interface, registerObject),
    interface_(interface.modelCast<dragModel, dispersedPhaseInterface>()),
    swarmCorrection_
    (
        dict.found("swarmCorrection")
      ? swarmCorrection::New(dict.subDict("swarmCorrection"), interface).ptr()
      : new swarmCorrections::noSwarm(dict, interface)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::dispersedDragModel::~dispersedDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::dispersedDragModel::Ki() const
{
    return
        0.75
       *CdRe()
       *swarmCorrection_->Cs()
       *interface_.continuous().rho()
       *interface_.continuous().thermo().nu()
       /sqr(interface_.dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::dispersedDragModel::K() const
{
    return
        max
        (
            interface_.dispersed(),
            interface_.dispersed().residualAlpha()
        )*Ki();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::dispersedDragModel::Kf() const
{
    return
        max
        (
            fvc::interpolate(interface_.dispersed()),
            interface_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki());
}


// ************************************************************************* //

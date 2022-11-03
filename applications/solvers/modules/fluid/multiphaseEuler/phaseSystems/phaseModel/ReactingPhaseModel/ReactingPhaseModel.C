/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "ReactingPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "combustionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ReactingPhaseModel<BasePhaseModel>::ReactingPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    reaction_(combustionModel::New(this->thermo_(), this->turbulence_()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ReactingPhaseModel<BasePhaseModel>::~ReactingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::ReactingPhaseModel<BasePhaseModel>::correctReactions()
{
    reaction_->correct();

    BasePhaseModel::correctReactions();
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix> Foam::ReactingPhaseModel<BasePhaseModel>::R
(
    volScalarField& Yi
) const
{
    return reaction_->R(Yi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::ReactingPhaseModel<BasePhaseModel>::Qdot() const
{
    return reaction_->Qdot();
}


// ************************************************************************* //

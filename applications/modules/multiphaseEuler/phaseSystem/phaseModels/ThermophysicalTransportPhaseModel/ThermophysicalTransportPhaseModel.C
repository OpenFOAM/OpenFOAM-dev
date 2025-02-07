/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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

#include "ThermophysicalTransportPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ThermophysicalTransportPhaseModel<BasePhaseModel>::
ThermophysicalTransportPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    thermophysicalTransport_
    (
        thermophysicalTransportModel::New
        (
            this->momentumTransport_,
            this->thermo_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ThermophysicalTransportPhaseModel<BasePhaseModel>::
~ThermophysicalTransportPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::ThermophysicalTransportPhaseModel<BasePhaseModel>::
predictThermophysicalTransport()
{
    BasePhaseModel::predictThermophysicalTransport();
    thermophysicalTransport_->predict();
}


template<class BasePhaseModel>
void Foam::ThermophysicalTransportPhaseModel<BasePhaseModel>::
correctThermophysicalTransport()
{
    BasePhaseModel::correctThermophysicalTransport();
    thermophysicalTransport_->correct();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::ThermophysicalTransportPhaseModel<BasePhaseModel>::kappaEff
(
    const label patchi
) const
{
    return thermophysicalTransport_->kappaEff(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermophysicalTransportPhaseModel<BasePhaseModel>::divq
(
    volScalarField& he
) const
{
    return thermophysicalTransport_->divq(he);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermophysicalTransportPhaseModel<BasePhaseModel>::divj
(
    volScalarField& Yi
) const
{
    return thermophysicalTransport_->divj(Yi);
}


// ************************************************************************* //

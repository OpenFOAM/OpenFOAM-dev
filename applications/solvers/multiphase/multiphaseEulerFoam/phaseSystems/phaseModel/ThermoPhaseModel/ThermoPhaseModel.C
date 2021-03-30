/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "ThermoPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoModel>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::ThermoPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    dynamicTransportModel(),
    thermo_(ThermoModel::New(fluid.mesh(), this->name()))
{
    thermo_->validate
    (
        IOobject::groupName(phaseModel::typeName, this->name()),
        "h",
        "e"
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoModel>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::~ThermoPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoModel>
bool Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::incompressible() const
{
    return thermo_().incompressible();
}


template<class BasePhaseModel, class ThermoModel>
bool Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::isochoric() const
{
    return thermo_().isochoric();
}


template<class BasePhaseModel, class ThermoModel>
const Foam::rhoThermo&
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::thermo() const
{
    return thermo_();
}


template<class BasePhaseModel, class ThermoModel>
Foam::rhoThermo&
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::thermoRef()
{
    return thermo_();
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::rho() const
{
    return thermo_->rho();
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::mu() const
{
    return thermo_->mu();
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::scalarField>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::mu
(
    const label patchi
) const
{
    return thermo_->mu(patchi);
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::nu() const
{
    return thermo_->nu();
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::scalarField>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::nu
(
    const label patchi
) const
{
    return thermo_->nu(patchi);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "SolidThermoPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoModel>
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::SolidThermoPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    viscosity(),
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
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::
~SolidThermoPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoModel>
bool
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::incompressible() const
{
    return thermo_().incompressible();
}


template<class BasePhaseModel, class ThermoModel>
bool Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::isochoric() const
{
    return thermo_().isochoric();
}


template<class BasePhaseModel, class ThermoModel>
const Foam::rhoThermo&
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::thermo() const
{
    return thermo_();
}


template<class BasePhaseModel, class ThermoModel>
Foam::rhoThermo&
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::thermo()
{
    return thermo_();
}


template<class BasePhaseModel, class ThermoModel>
const Foam::rhoFluidThermo&
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::fluidThermo() const
{
    NotImplemented;
    return refCast<const rhoFluidThermo>(thermo_());
}


template<class BasePhaseModel, class ThermoModel>
Foam::rhoFluidThermo&
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::fluidThermo()
{
    NotImplemented;
    return refCast<rhoFluidThermo>(thermo_());
}


template<class BasePhaseModel, class ThermoModel>
const Foam::volScalarField&
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::rho() const
{
    return thermo_->rho();
}


template<class BasePhaseModel, class ThermoModel>
Foam::volScalarField&
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::rho()
{
    return thermo_->rho();
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::mu() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::mu
(
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::nu() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoPhaseModel<BasePhaseModel, ThermoModel>::nu
(
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


// ************************************************************************* //

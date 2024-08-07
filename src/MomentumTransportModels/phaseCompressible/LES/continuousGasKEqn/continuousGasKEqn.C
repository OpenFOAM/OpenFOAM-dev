/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

#include "continuousGasKEqn.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
continuousGasKEqn<BasicMomentumTransportModel>::continuousGasKEqn
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    kEqn<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),

    liquidTurbulencePtr_(nullptr),

    alphaInversion_("alphaInversion", this->coeffDict(), 0.7)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool continuousGasKEqn<BasicMomentumTransportModel>::read()
{
    if (kEqn<BasicMomentumTransportModel>::read())
    {
        alphaInversion_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
const momentumTransportModel&
continuousGasKEqn<BasicMomentumTransportModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const phaseModel& gas = refCast<const phaseModel>(this->properties());
        const phaseSystem& fluid = gas.fluid();
        const phaseModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
            &U.db().lookupType<momentumTransportModel>(liquid.name());
    }

    return *liquidTurbulencePtr_;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
continuousGasKEqn<BasicMomentumTransportModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const momentumTransportModel& liquidTurbulence = this->liquidTurbulence();

    return
    (
        max(alphaInversion_ - alpha, scalar(0))
       *rho
       *min
        (
            this->Ce_*sqrt(liquidTurbulence.k())/this->delta(),
            1.0/U.time().deltaT()
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
continuousGasKEqn<BasicMomentumTransportModel>::kSource() const
{
    const momentumTransportModel& liquidTurbulence = this->liquidTurbulence();
    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        phaseTransferCoeff*liquidTurbulence.k()
      - fvm::Sp(phaseTransferCoeff, this->k_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //

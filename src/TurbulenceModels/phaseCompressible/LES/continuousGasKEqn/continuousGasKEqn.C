/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
continuousGasKEqn<BasicTurbulenceModel>::continuousGasKEqn
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEqn<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    liquidTurbulencePtr_(NULL),

    alphaInversion_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaInversion",
            this->coeffDict_,
            0.7
        )
    )
{
    if (type == typeName)
    {
        kEqn<BasicTurbulenceModel>::correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool continuousGasKEqn<BasicTurbulenceModel>::read()
{
    if (kEqn<BasicTurbulenceModel>::read())
    {
        alphaInversion_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
const turbulenceModel&
continuousGasKEqn<BasicTurbulenceModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& gas = this->transport();
        const twoPhaseSystem& fluid = gas.fluid();
        const transportModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &U.db().lookupObject<turbulenceModel>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    liquid.name()
                )
            );
    }

    return *liquidTurbulencePtr_;
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
continuousGasKEqn<BasicTurbulenceModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const turbulenceModel& liquidTurbulence = this->liquidTurbulence();

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


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
continuousGasKEqn<BasicTurbulenceModel>::kSource() const
{
    const turbulenceModel& liquidTurbulence = this->liquidTurbulence();
    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        phaseTransferCoeff*liquidTurbulence.k()
      - fvm::Sp(phaseTransferCoeff, this->k_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //

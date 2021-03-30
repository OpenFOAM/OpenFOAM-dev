/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "NicenoKEqn.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "phaseSystem.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
NicenoKEqn<BasicMomentumTransportModel>::NicenoKEqn
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
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
        transport,
        type
    ),

    gasTurbulencePtr_(nullptr),

    alphaInversion_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaInversion",
            this->coeffDict_,
            0.3
        )
    ),

    Cp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp",
            this->coeffDict_,
            this->Ck_.value()
        )
    ),

    Cmub_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmub",
            this->coeffDict_,
            0.6
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool NicenoKEqn<BasicMomentumTransportModel>::read()
{
    if (kEqn<BasicMomentumTransportModel>::read())
    {
        alphaInversion_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());
        Cmub_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
const PhaseCompressibleMomentumTransportModel
<
    typename BasicMomentumTransportModel::transportModel
>&
NicenoKEqn<BasicMomentumTransportModel>::gasTurbulence() const
{
    if (!gasTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const phaseModel& liquid = refCast<const phaseModel>(this->transport());
        const phaseSystem& fluid = liquid.fluid();
        const phaseModel& gas = fluid.otherPhase(liquid);

        gasTurbulencePtr_ =
           &U.db().lookupObject
            <
                PhaseCompressibleMomentumTransportModel<transportModel>
            >
            (
                IOobject::groupName
                (
                    momentumTransportModel::typeName,
                    gas.name()
                )
            );
    }

    return *gasTurbulencePtr_;
}


template<class BasicMomentumTransportModel>
void NicenoKEqn<BasicMomentumTransportModel>::correctNut()
{
    const PhaseCompressibleMomentumTransportModel<transportModel>&
        gasTurbulence = this->gasTurbulence();

    const phaseModel& liquid = refCast<const phaseModel>(this->transport());
    const phaseSystem& fluid = liquid.fluid();
    const phaseModel& gas = fluid.otherPhase(liquid);

    this->nut_ =
        this->Ck_*sqrt(this->k_)*this->delta()
      + Cmub_*gas.d()*gasTurbulence.alpha()
       *(mag(this->U_ - gasTurbulence.U()));

    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> NicenoKEqn<BasicMomentumTransportModel>::bubbleG() const
{
    const PhaseCompressibleMomentumTransportModel<transportModel>&
        gasTurbulence = this->gasTurbulence();

    const phaseModel& liquid = refCast<const phaseModel>(this->transport());
    const phaseSystem& fluid = liquid.fluid();
    const phaseModel& gas = fluid.otherPhase(liquid);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(this->U_ - gasTurbulence.U()));

    tmp<volScalarField> bubbleG
    (
        Cp_*sqr(magUr)*drag.K()/liquid.rho()
    );

    return bubbleG;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
NicenoKEqn<BasicMomentumTransportModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const momentumTransportModel& gasTurbulence = this->gasTurbulence();

    return
    (
        max(alphaInversion_ - alpha, scalar(0))
       *rho
       *min
        (
            this->Ce_*sqrt(gasTurbulence.k())/this->delta(),
            1.0/U.time().deltaT()
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> NicenoKEqn<BasicMomentumTransportModel>::kSource() const
{
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const PhaseCompressibleMomentumTransportModel<transportModel>&
        gasTurbulence = this->gasTurbulence();

    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        alpha*rho*bubbleG()
      + phaseTransferCoeff*gasTurbulence.k()
      - fvm::Sp(phaseTransferCoeff, this->k_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //

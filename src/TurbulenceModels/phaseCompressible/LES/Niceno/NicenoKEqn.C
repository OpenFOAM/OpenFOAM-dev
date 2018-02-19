/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
#include "fvOptions.H"
#include "twoPhaseSystem.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
NicenoKEqn<BasicTurbulenceModel>::NicenoKEqn
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

template<class BasicTurbulenceModel>
bool NicenoKEqn<BasicTurbulenceModel>::read()
{
    if (kEqn<BasicTurbulenceModel>::read())
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


template<class BasicTurbulenceModel>
const PhaseCompressibleTurbulenceModel
<
    typename BasicTurbulenceModel::transportModel
>&
NicenoKEqn<BasicTurbulenceModel>::gasTurbulence() const
{
    if (!gasTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& liquid = this->transport();
        const twoPhaseSystem& fluid =
            refCast<const twoPhaseSystem>(liquid.fluid());
        const transportModel& gas = fluid.otherPhase(liquid);

        gasTurbulencePtr_ =
           &U.db()
           .lookupObject<PhaseCompressibleTurbulenceModel<transportModel>>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    gas.name()
                )
            );
    }

    return *gasTurbulencePtr_;
}


template<class BasicTurbulenceModel>
void NicenoKEqn<BasicTurbulenceModel>::correctNut()
{
    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    this->nut_ =
        this->Ck_*sqrt(this->k_)*this->delta()
      + Cmub_*gasTurbulence.transport().d()*gasTurbulence.alpha()
       *(mag(this->U_ - gasTurbulence.U()));

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> NicenoKEqn<BasicTurbulenceModel>::bubbleG() const
{
    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    const transportModel& liquid = this->transport();
    const twoPhaseSystem& fluid =
        refCast<const twoPhaseSystem>(liquid.fluid());
    const transportModel& gas = fluid.otherPhase(liquid);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(this->U_ - gasTurbulence.U()));

    tmp<volScalarField> bubbleG
    (
        Cp_*sqr(magUr)*drag.K()/liquid.rho()
    );

    return bubbleG;
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
NicenoKEqn<BasicTurbulenceModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const turbulenceModel& gasTurbulence = this->gasTurbulence();

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


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> NicenoKEqn<BasicTurbulenceModel>::kSource() const
{
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

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

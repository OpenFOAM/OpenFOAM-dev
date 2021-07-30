/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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

#include "kOmegaSSTSato.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "phaseSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmegaSSTSato<BasicMomentumTransportModel>::kOmegaSSTSato
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
    kOmegaSST<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),

    phase_(refCast<const phaseModel>(viscosity)),

    hasDispersedPhaseNames_(this->coeffDict_.found("dispersedPhases")),

    dispersedPhaseNames_
    (
        this->coeffDict_.lookupOrDefault("dispersedPhases", hashedWordList())
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
bool kOmegaSSTSato<BasicMomentumTransportModel>::read()
{
    if (kOmegaSST<BasicMomentumTransportModel>::read())
    {
        Cmub_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
UPtrList<const phaseModel>
kOmegaSSTSato<BasicMomentumTransportModel>::getDispersedPhases() const
{
    UPtrList<const phaseModel> dispersedPhases;

    const phaseSystem& fluid = phase_.fluid();

    if (hasDispersedPhaseNames_)
    {
        dispersedPhases.resize(dispersedPhaseNames_.size());

        forAll(dispersedPhaseNames_, dispersedPhasei)
        {
            dispersedPhases.set
            (
                dispersedPhasei,
                &fluid.phases()[dispersedPhaseNames_[dispersedPhasei]]
            );
        }
    }
    else
    {
        dispersedPhases.resize(fluid.movingPhases().size() - 1);

        label dispersedPhasei = 0;

        forAll(fluid.movingPhases(), movingPhasei)
        {
            const phaseModel& otherPhase = fluid.movingPhases()[movingPhasei];

            if (&otherPhase != &phase_)
            {
                dispersedPhases.set
                (
                    dispersedPhasei ++,
                    &otherPhase
                );
            }
        }
    }

    return dispersedPhases;
}


template<class BasicMomentumTransportModel>
void kOmegaSSTSato<BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    volScalarField yPlus
    (
        pow(this->betaStar_, 0.25)*this->y_*sqrt(this->k_)/this->nu()
    );

    this->nut_ =
        this->a1_*this->k_/max(this->a1_*this->omega_, this->b1_*F2*sqrt(S2));

    UPtrList<const phaseModel> dispersedPhases(getDispersedPhases());
    forAllIter(UPtrList<const phaseModel>, dispersedPhases, phaseIter)
    {
        this->nut_ +=
            sqr(1 - exp(-yPlus/16.0))
           *Cmub_*phaseIter().d()*phaseIter()
           *(mag(this->U_ - phaseIter().U()));
    }

    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void kOmegaSSTSato<BasicMomentumTransportModel>::correct()
{
    kOmegaSST<BasicMomentumTransportModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

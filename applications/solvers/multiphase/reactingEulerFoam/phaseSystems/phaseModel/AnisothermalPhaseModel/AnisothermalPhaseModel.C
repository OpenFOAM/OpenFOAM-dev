/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "AnisothermalPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::AnisothermalPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    K_
    (
        IOobject
        (
            IOobject::groupName("K", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("K", sqr(dimVelocity), scalar(0))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::~AnisothermalPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::AnisothermalPhaseModel<BasePhaseModel>::compressible() const
{
    return !this->thermo().incompressible();
}


template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::correctKinematics()
{
    BasePhaseModel::correctKinematics();
    K_ = 0.5*magSqr(this->U());
}


template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    this->thermo_->correct();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::filterPressureWork
(
    const tmp<volScalarField>& pressureWork
) const
{
    const volScalarField& alpha = *this;

    scalar pressureWorkAlphaLimit =
        this->thermo_->lookupOrDefault("pressureWorkAlphaLimit", 0.0);

    if (pressureWorkAlphaLimit > 0)
    {
        return
        (
            max(alpha - pressureWorkAlphaLimit, scalar(0))
           /max(alpha - pressureWorkAlphaLimit, pressureWorkAlphaLimit)
        )*pressureWork;
    }
    else
    {
        return pressureWork;
    }
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;
    const volVectorField& U = this->U();
    const surfaceScalarField& alphaPhi = this->alphaPhi();
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi();

    const volScalarField& contErr(this->continuityError());

    const volScalarField alphaEff(this->alphaEff());

    volScalarField& he = this->thermo_->he();

    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha, this->rho(), he)
      + fvm::div(alphaRhoPhi, he)
      - fvm::Sp(contErr, he)

      + fvc::ddt(alpha, this->rho(), K_) + fvc::div(alphaRhoPhi, K_)
      - contErr*K_

      - fvm::laplacian
        (
            fvc::interpolate(alpha)
           *fvc::interpolate(alphaEff),
            he
        )
     ==
        alpha*this->Qdot()
    );

    // Add the appropriate pressure-work term
    if (he.name() == this->thermo_->phasePropertyName("e"))
    {
        tEEqn.ref() += filterPressureWork
        (
            fvc::div(fvc::absolute(alphaPhi, alpha, U), this->thermo().p())
          + this->thermo().p()*fvc::ddt(alpha)
        );
    }
    else if (this->thermo_->dpdt())
    {
        tEEqn.ref() -= filterPressureWork(alpha*this->fluid().dpdt());
    }

    return tEEqn;
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::AnisothermalPhaseModel<BasePhaseModel>::K() const
{
    return K_;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseName),
    D_
    (
        IOobject
        (
            IOobject::groupName("D", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("D", dimless/dimTime, 0)
    ),
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
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;
    const surfaceScalarField& alphaPhi = this->alphaPhi();
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi();

    const volScalarField& contErr(this->continuityError());

    const volScalarField alphaEff
    (
        this->thermo_->alphaEff(this->turbulence().mut())
    );

    volScalarField& he = this->thermo_->he();

    return
    (
        fvm::ddt(alpha, this->rho(), he) + fvm::div(alphaRhoPhi, he)
      - fvm::Sp(contErr, he)

      + fvc::ddt(alpha, this->rho(), K_) + fvc::div(alphaRhoPhi, K_)
      - contErr*K_

      - fvm::laplacian
        (
            fvc::interpolate(alpha)
           *fvc::interpolate(alphaEff),
            he
        )

      + (
            he.name() == this->thermo_->phasePropertyName("e")
          ? fvc::ddt(alpha)*this->thermo().p()
          + fvc::div(alphaPhi, this->thermo().p())
          : -alpha*this->fluid().dpdt()
        )
     ==
        this->Sh()
    );
}


template<class BasePhaseModel>
bool Foam::AnisothermalPhaseModel<BasePhaseModel>::compressible() const
{
    return true;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::D() const
{
    return D_;
}


template<class BasePhaseModel>
void
Foam::AnisothermalPhaseModel<BasePhaseModel>::D(const volScalarField& D)
{
    D_ = D;
}


// ************************************************************************* //

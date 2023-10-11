/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::filterPressureWork
(
    const tmp<volScalarField>& pressureWork
) const
{
    const volScalarField& alpha = *this;

    scalar pressureWorkAlphaLimit =
        this->thermo_->properties()
       .lookupOrDefault("pressureWorkAlphaLimit", 0.0);

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::AnisothermalPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    g_(fluid.mesh().lookupObject<uniformDimensionedVectorField>("g")),
    thermophysicalTransport_
    (
        PhaseThermophysicalTransportModel
        <
            phaseCompressible::momentumTransportModel,
            transportThermoModel
        >::New(this->momentumTransport_, this->thermo_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::~AnisothermalPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    this->thermo_->correct();
}


template<class BasePhaseModel>
bool Foam::AnisothermalPhaseModel<BasePhaseModel>::isothermal() const
{
    return false;
}


template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::
predictThermophysicalTransport()
{
    BasePhaseModel::predictThermophysicalTransport();
    thermophysicalTransport_->predict();
}


template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::
correctThermophysicalTransport()
{
    BasePhaseModel::correctThermophysicalTransport();
    thermophysicalTransport_->correct();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::kappaEff(const label patchi) const
{
    return thermophysicalTransport_->kappaEff(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::divq(volScalarField& he) const
{
    return thermophysicalTransport_->divq(he);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::divj(volScalarField& Yi) const
{
    return thermophysicalTransport_->divj(Yi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;
    const volScalarField& rho = this->rho();

    const tmp<volVectorField> tU(this->U());
    const volVectorField& U(tU());

    const tmp<surfaceScalarField> talphaRhoPhi(this->alphaRhoPhi());
    const surfaceScalarField& alphaRhoPhi(talphaRhoPhi());

    const tmp<volScalarField> tcontErr(this->continuityError());
    const volScalarField& contErr(tcontErr());

    tmp<volScalarField> tK(this->K());
    const volScalarField& K(tK());

    volScalarField& he = this->thermo_->he();

    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha, rho, he)
      + fvm::div(alphaRhoPhi, he)
      - fvm::Sp(contErr, he)

      + fvc::ddt(alpha, rho, K) + fvc::div(alphaRhoPhi, K)
      - contErr*K

      + this->divq(he)
     ==
        alpha*rho*(U&g_)
      + alpha*this->Qdot()
    );

    // Add the appropriate pressure-work term
    if (he.name() == this->thermo_->phasePropertyName("e"))
    {
        tEEqn.ref() += filterPressureWork
        (
            fvc::div
            (
                fvc::absolute(alphaRhoPhi, alpha, rho, U),
                this->fluidThermo().p()/rho
            )
          + (fvc::ddt(alpha) - contErr/rho)*this->fluidThermo().p()
        );
    }
    else if (this->thermo_->dpdt())
    {
        tEEqn.ref() -= filterPressureWork(alpha*this->fluid().dpdt());
    }

    return tEEqn;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "compressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleInterPhaseTransportModel::compressibleInterPhaseTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& rhoPhi,
    const surfaceScalarField& alphaPhi1,
    const compressibleTwoPhaseMixture& mixture
)
:
    twoPhaseTransport_(false),
    mixture_(mixture),
    phi_(phi),
    alphaPhi1_(alphaPhi1)
{
    {
        IOdictionary momentumTransport
        (
            IOobject
            (
                momentumTransportModel::typeName,
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word simulationType
        (
            momentumTransport.lookup("simulationType")
        );

        if (simulationType == "twoPhaseTransport")
        {
            twoPhaseTransport_ = true;
        }
    }

    if (twoPhaseTransport_)
    {
        const volScalarField& alpha1(mixture_.alpha1());
        const volScalarField& alpha2(mixture_.alpha2());

        const volScalarField& rho1 = mixture_.thermo1().rho();
        const volScalarField& rho2 = mixture_.thermo2().rho();

        alphaRhoPhi1_ =
        (
            new surfaceScalarField
            (
                IOobject::groupName("alphaRhoPhi", alpha1.group()),
                fvc::interpolate(rho1)*alphaPhi1_
            )
        );

        alphaRhoPhi2_ =
        (
            new surfaceScalarField
            (
                IOobject::groupName("alphaRhoPhi", alpha2.group()),
                fvc::interpolate(rho2)*(phi_ - alphaPhi1_)
            )
        );

        turbulence1_ =
        (
            phaseCompressible::momentumTransportModel::New
            (
                alpha1,
                rho1,
                U,
                alphaRhoPhi1_(),
                phi,
                mixture.thermo1()
            )
        );

        turbulence2_ =
        (
            phaseCompressible::momentumTransportModel::New
            (
                alpha2,
                rho2,
                U,
                alphaRhoPhi2_(),
                phi,
                mixture.thermo2()
            )
        );
    }
    else
    {
        turbulence_ = compressible::momentumTransportModel::New
        (
            rho,
            U,
            rhoPhi,
            mixture
        );

        turbulence_->validate();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseTransportModel::alphaEff() const
{
    /* ***HGW
    if (twoPhaseTransport_)
    {
        return
            mixture_.alpha1()*mixture_.thermo1().alphaEff
            (
                turbulence1_->alphat()
            )
          + mixture_.alpha2()*mixture_.thermo2().alphaEff
            (
                turbulence2_->alphat()
            );
    }
    else
    {
        return mixture_.alphaEff(turbulence_->alphat());
    }
    */

    if (twoPhaseTransport_)
    {
        return
            mixture_.alpha1()*mixture_.thermo1().alphaEff
            (
                mixture_.thermo1().rho()*turbulence1_->nut()
            )
          + mixture_.alpha2()*mixture_.thermo2().alphaEff
            (
                mixture_.thermo2().rho()*turbulence2_->nut()
            );
    }
    else
    {
        const volScalarField alphat(mixture_.rho()*turbulence_->nut());

        return
            mixture_.alpha1()*mixture_.thermo1().alphaEff(alphat)
          + mixture_.alpha2()*mixture_.thermo2().alphaEff(alphat);
    }
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::compressibleInterPhaseTransportModel::divDevTau
(
    volVectorField& U
) const
{
    if (twoPhaseTransport_)
    {
        return
            turbulence1_->divDevTau(U)
          + turbulence2_->divDevTau(U);
    }
    else
    {
        return turbulence_->divDevTau(U);
    }
}


void Foam::compressibleInterPhaseTransportModel::correctPhasePhi()
{
    if (twoPhaseTransport_)
    {
        const volScalarField& rho1 = mixture_.thermo1().rho();
        const volScalarField& rho2 = mixture_.thermo2().rho();

        alphaRhoPhi1_.ref() = fvc::interpolate(rho1)*alphaPhi1_;
        alphaRhoPhi2_.ref() = fvc::interpolate(rho2)*(phi_ - alphaPhi1_);
    }
}


void Foam::compressibleInterPhaseTransportModel::correct()
{
    if (twoPhaseTransport_)
    {
        turbulence1_->correct();
        turbulence2_->correct();
    }
    else
    {
        turbulence_->correct();
    }
}


// ************************************************************************* //

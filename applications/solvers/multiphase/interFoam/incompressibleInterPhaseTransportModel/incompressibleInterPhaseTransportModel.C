/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "incompressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleInterPhaseTransportModel::
incompressibleInterPhaseTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& alphaPhi10,
    const incompressibleTwoPhaseMixture& mixture
)
:
    twoPhaseTransport_(false),
    mixture_(mixture),
    phi_(phi),
    alphaPhi10_(alphaPhi10)
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

        alphaPhi2_ =
        (
            new surfaceScalarField
            (
                IOobject::groupName("alphaPhi", alpha2.group()),
                (phi_ - alphaPhi10_)
            )
        );

        turbulence1_ =
        (
            phaseIncompressible::momentumTransportModel::New
            (
                alpha1,
                U,
                alphaPhi10_,
                phi,
                mixture.nuModel1()
            )
        );

        turbulence2_ =
        (
            phaseIncompressible::momentumTransportModel::New
            (
                alpha2,
                U,
                alphaPhi2_(),
                phi,
                mixture.nuModel2()
            )
        );
    }
    else
    {
        turbulence_ = incompressible::momentumTransportModel::New
        (
            U,
            phi,
            mixture
        );

        turbulence_->validate();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleInterPhaseTransportModel::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    if (twoPhaseTransport_)
    {
        return
          mixture_.rho1()*turbulence1_->divDevSigma(U)
        + mixture_.rho2()*turbulence2_->divDevSigma(U);
    }
    else
    {
        return turbulence_->divDevTau(rho, U);
    }
}


void Foam::incompressibleInterPhaseTransportModel::correctPhasePhi()
{
    if (twoPhaseTransport_)
    {
        alphaPhi2_.ref() = (phi_ - alphaPhi10_);
    }
}


void Foam::incompressibleInterPhaseTransportModel::correct()
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

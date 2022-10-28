/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "compressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixture::compressibleTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    twoPhaseMixture(U.mesh()),
    interfaceProperties(alpha1(), alpha2(), U, *this),

    totalInternalEnergy_
    (
        lookupOrDefault<Switch>("totalInternalEnergy", true)
    ),

    p_
    (
        IOobject
        (
            "p",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    T_
    (
        IOobject
        (
            "T",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    thermo1_(nullptr),
    thermo2_(nullptr),
    rho_
    (
        IOobject
        (
            "rho",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("rho", dimDensity, 0)
    ),
    Alpha1_
    (
        IOobject
        (
            IOobject::groupName("Alpha", phase1Name()),
            U.mesh().time().timeName(),
            U.mesh()
        ),
        alpha1(),
        calculatedFvPatchScalarField::typeName
    ),
    Alpha2_
    (
        IOobject
        (
            IOobject::groupName("Alpha", phase2Name()),
            U.mesh().time().timeName(),
            U.mesh()
        ),
        alpha2(),
        calculatedFvPatchScalarField::typeName
    )
{
    {
        volScalarField T1
        (
            IOobject
            (
                IOobject::groupName("T", phase1Name()),
                U.mesh().time().timeName(),
                U.mesh()
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );
        T1.write();
    }

    {
        volScalarField T2
        (
            IOobject
            (
                IOobject::groupName("T", phase2Name()),
                U.mesh().time().timeName(),
                U.mesh()
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );
        T2.write();
    }

    // Note: we're writing files to be read in immediately afterwards.
    //       Avoid any thread-writing problems.
    // fileHandler().flush();

    thermo1_ = rhoThermo::New(U.mesh(), phase1Name());
    thermo2_ = rhoThermo::New(U.mesh(), phase2Name());

    // thermo1_->validate(phase1Name(), "e");
    // thermo2_->validate(phase2Name(), "e");

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixture::~compressibleTwoPhaseMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibleTwoPhaseMixture::correctThermo()
{
    thermo1_->T() = T_;
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->T() = T_;
    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();
}


void Foam::compressibleTwoPhaseMixture::correct()
{
    const volScalarField alphaRho1(alpha1()*thermo1_->rho());
    const volScalarField alphaRho2(alpha2()*thermo2_->rho());

    rho_ = alphaRho1 + alphaRho2;
    Alpha1_ = alphaRho1/rho_;
    Alpha2_ = alphaRho2/rho_;

    interfaceProperties::correct();
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixture::nu() const
{
    return (alpha1()*thermo1_->mu() + alpha2()*thermo2_->mu())/rho_;
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixture::nu
(
    const label patchi
) const
{
    return
    (
        alpha1().boundaryField()[patchi]*thermo1_->mu(patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->mu(patchi)
    )/rho_.boundaryField()[patchi];
}


bool Foam::compressibleTwoPhaseMixture::read()
{
    if (twoPhaseMixture::read())
    {
        totalInternalEnergy_ =
            lookupOrDefault<Switch>("totalInternalEnergy", true);

        return interfaceProperties::read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

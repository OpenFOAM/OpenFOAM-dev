/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "SchnerrSauer.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseChangeModels
{
    defineTypeNameAndDebug(SchnerrSauer, 0);
    addToRunTimeSelectionTable
    (
        cavitationModel,
        SchnerrSauer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::SchnerrSauer::SchnerrSauer
(
    const twoPhaseMixtureThermo& mixture
)
:
    cavitationModel(typeName, mixture),

    n_("n", dimless/dimVolume, twoPhaseChangeModelCoeffs_),
    dNuc_("dNuc", dimLength, twoPhaseChangeModelCoeffs_),
    Cc_("Cc", dimless, twoPhaseChangeModelCoeffs_),
    Cv_("Cv", dimless, twoPhaseChangeModelCoeffs_),

    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::twoPhaseChangeModels::SchnerrSauer::rRb
(
    const volScalarField::Internal& limitedAlpha1
) const
{
    return pow
    (
        ((4*constant::mathematical::pi*n_)/3)
       *limitedAlpha1/(1.0 + alphaNuc() - limitedAlpha1),
        1.0/3.0
    );
}


Foam::dimensionedScalar
Foam::twoPhaseChangeModels::SchnerrSauer::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::twoPhaseChangeModels::SchnerrSauer::pCoeff
(
    const volScalarField::Internal& p
) const
{
    const volScalarField::Internal limitedAlpha1
    (
        min(max(mixture_.alpha1()(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal rho
    (
        limitedAlpha1*rho1()
      + (scalar(1) - limitedAlpha1)*rho2()
    );

    return
        (3*rho1()*rho2())*sqrt(2/(3*rho1()))
       *rRb(limitedAlpha1)/(rho*sqrt(mag(p - pSat()) + 0.01*pSat()));
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::SchnerrSauer::mDotAlphal() const
{
    const volScalarField::Internal& p =
        mixture_.alpha1().db().lookupObject<volScalarField>("p");

    const volScalarField::Internal pCoeff(this->pCoeff(p));

    const volScalarField::Internal limitedAlpha1
    (
        min(max(mixture_.alpha1()(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField::Internal>>
    (
        Cc_*limitedAlpha1*pCoeff*max(p - pSat(), p0_),

        Cv_*(1.0 + alphaNuc() - limitedAlpha1)*pCoeff*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::SchnerrSauer::mDotP() const
{
    const volScalarField::Internal& p =
        mixture_.alpha1().db().lookupObject<volScalarField>("p");

    const volScalarField::Internal pCoeff(this->pCoeff(p));

    const volScalarField::Internal limitedAlpha1
    (
        min(max(mixture_.alpha1()(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal apCoeff(limitedAlpha1*pCoeff);

    return Pair<tmp<volScalarField::Internal>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos0(p - pSat())*apCoeff,

        (-Cv_)*(1.0 + alphaNuc() - limitedAlpha1)*neg(p - pSat())*apCoeff
    );
}


void Foam::twoPhaseChangeModels::SchnerrSauer::correct()
{
    cavitationModel::correct();
}


bool Foam::twoPhaseChangeModels::SchnerrSauer::read()
{
    if (cavitationModel::read())
    {
        twoPhaseChangeModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        twoPhaseChangeModelCoeffs_.lookup("n") >> n_;
        twoPhaseChangeModelCoeffs_.lookup("dNuc") >> dNuc_;
        twoPhaseChangeModelCoeffs_.lookup("Cc") >> Cc_;
        twoPhaseChangeModelCoeffs_.lookup("Cv") >> Cv_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

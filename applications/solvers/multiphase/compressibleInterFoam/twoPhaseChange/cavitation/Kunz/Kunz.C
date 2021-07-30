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

#include "Kunz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseChangeModels
{
    defineTypeNameAndDebug(Kunz, 0);
    addToRunTimeSelectionTable(cavitationModel, Kunz, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::Kunz::Kunz
(
    const compressibleTwoPhaseMixture& mixture
)
:
    cavitationModel(typeName, mixture),

    UInf_("UInf", dimVelocity, twoPhaseChangeModelCoeffs_),
    tInf_("tInf", dimTime, twoPhaseChangeModelCoeffs_),
    Cc_("Cc", dimless, twoPhaseChangeModelCoeffs_),
    Cv_("Cv", dimless, twoPhaseChangeModelCoeffs_),

    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::Kunz::mDotAlphal() const
{
    const volScalarField::Internal& p =
        mixture_.alpha1().db().lookupObject<volScalarField>("p");

    const volScalarField::Internal mcCoeff_(Cc_*rho2()/tInf_);
    const volScalarField::Internal mvCoeff_
    (
        Cv_*rho2()/(0.5*rho1()*sqr(UInf_)*tInf_)
    );

    const volScalarField::Internal limitedAlpha1
    (
        min(max(mixture_.alpha1()(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff_*sqr(limitedAlpha1)
       *max(p - pSat(), p0_)/max(p - pSat(), 0.01*pSat()),

        mvCoeff_*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::Kunz::mDotP() const
{
    const volScalarField::Internal& p =
        mixture_.alpha1().db().lookupObject<volScalarField>("p");

    const volScalarField::Internal mcCoeff_(Cc_*rho2()/tInf_);
    const volScalarField::Internal mvCoeff_
    (
        Cv_*rho2()/(0.5*rho1()*sqr(UInf_)*tInf_)
    );

    const volScalarField::Internal limitedAlpha1
    (
        min(max(mixture_.alpha1()(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff_*sqr(limitedAlpha1)*(1.0 - limitedAlpha1)
       *pos0(p - pSat())/max(p - pSat(), 0.01*pSat()),

        (-mvCoeff_)*limitedAlpha1*neg(p - pSat())
    );
}


void Foam::twoPhaseChangeModels::Kunz::correct()
{
    cavitationModel::correct();
}


bool Foam::twoPhaseChangeModels::Kunz::read()
{
    if (cavitationModel::read())
    {
        twoPhaseChangeModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        twoPhaseChangeModelCoeffs_.lookup("UInf") >> UInf_;
        twoPhaseChangeModelCoeffs_.lookup("tInf") >> tInf_;
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

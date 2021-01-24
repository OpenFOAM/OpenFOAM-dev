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

#include "Merkle.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseChangeModels
{
    defineTypeNameAndDebug(Merkle, 0);
    addToRunTimeSelectionTable(cavitationModel, Merkle, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::Merkle::Merkle
(
    const immiscibleIncompressibleTwoPhaseMixture& mixture
)
:
    cavitationModel(typeName, mixture),

    UInf_("UInf", dimVelocity, twoPhaseChangeModelCoeffs_),
    tInf_("tInf", dimTime, twoPhaseChangeModelCoeffs_),
    Cc_("Cc", dimless, twoPhaseChangeModelCoeffs_),
    Cv_("Cv", dimless, twoPhaseChangeModelCoeffs_),

    p0_("0", pSat().dimensions(), 0.0),

    mcCoeff_(Cc_/(0.5*sqr(UInf_)*tInf_)),
    mvCoeff_(Cv_*mixture_.rho1()/(0.5*sqr(UInf_)*tInf_*mixture_.rho2()))
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::Merkle::mDotAlphal() const
{
    const volScalarField& p =
        mixture_.U().db().lookupObject<volScalarField>("p");

    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*max(p - pSat(), p0_),
        mvCoeff_*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::Merkle::mDotP() const
{
    const volScalarField& p =
        mixture_.U().db().lookupObject<volScalarField>("p");

    const volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*(1.0 - limitedAlpha1)*pos0(p - pSat()),
        (-mvCoeff_)*limitedAlpha1*neg(p - pSat())
    );
}


void Foam::twoPhaseChangeModels::Merkle::correct()
{
    cavitationModel::correct();
}


bool Foam::twoPhaseChangeModels::Merkle::read()
{
    if (cavitationModel::read())
    {
        twoPhaseChangeModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        twoPhaseChangeModelCoeffs_.lookup("UInf") >> UInf_;
        twoPhaseChangeModelCoeffs_.lookup("tInf") >> tInf_;
        twoPhaseChangeModelCoeffs_.lookup("Cc") >> Cc_;
        twoPhaseChangeModelCoeffs_.lookup("Cv") >> Cv_;

        mcCoeff_ = Cc_/(0.5*sqr(UInf_)*tInf_);
        mvCoeff_ = Cv_*mixture_.rho1()/(0.5*sqr(UInf_)*tInf_*mixture_.rho2());

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

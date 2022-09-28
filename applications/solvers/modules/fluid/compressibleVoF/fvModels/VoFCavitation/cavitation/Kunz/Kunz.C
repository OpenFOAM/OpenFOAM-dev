/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
namespace compressible
{
namespace cavitationModels
{
    defineTypeNameAndDebug(Kunz, 0);
    addToRunTimeSelectionTable(cavitationModel, Kunz, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressible::cavitationModels::Kunz::Kunz
(
    const dictionary& dict,
    const compressibleTwoPhaseMixture& mixture
)
:
    cavitationModel(dict, mixture),

    UInf_("UInf", dimVelocity, dict),
    tInf_("tInf", dimTime, dict),
    Cc_("Cc", dimless, dict),
    Cv_("Cv", dimless, dict),

    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::compressible::cavitationModels::Kunz::mDotAlphal() const
{
    const volScalarField::Internal& p =
        mixture_.alpha1().db().lookupObject<volScalarField>("p");

    const volScalarField::Internal mcCoeff_(Cc_*mixture_.rho2()()/tInf_);
    const volScalarField::Internal mvCoeff_
    (
        Cv_*mixture_.rho2()()/(0.5*mixture_.rho1()()*sqr(UInf_)*tInf_)
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
Foam::compressible::cavitationModels::Kunz::mDotP() const
{
    const volScalarField::Internal& p =
        mixture_.alpha1().db().lookupObject<volScalarField>("p");

    const volScalarField::Internal mcCoeff_(Cc_*mixture_.rho2()()/tInf_);
    const volScalarField::Internal mvCoeff_
    (
        Cv_*mixture_.rho2()()/(0.5*mixture_.rho1()()*sqr(UInf_)*tInf_)
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


void Foam::compressible::cavitationModels::Kunz::correct()
{}


bool Foam::compressible::cavitationModels::Kunz::read
(
    const dictionary& dict
)
{
    if (cavitationModel::read(dict))
    {
        dict.lookup("UInf") >> UInf_;
        dict.lookup("tInf") >> tInf_;
        dict.lookup("Cc") >> Cc_;
        dict.lookup("Cv") >> Cv_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

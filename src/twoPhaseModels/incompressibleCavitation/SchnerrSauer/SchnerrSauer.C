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

#include "SchnerrSauer.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cavitationModels
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

Foam::cavitationModels::SchnerrSauer::SchnerrSauer
(
    const dictionary& dict,
    const incompressibleTwoPhases& phases
)
:
    cavitationModel(dict, phases),

    n_("n", dimless/dimVolume, dict),
    dNuc_("dNuc", dimLength, dict),
    Cc_("Cc", dimless, dict),
    Cv_("Cv", dimless, dict),

    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::cavitationModels::SchnerrSauer::rRb
(
    const volScalarField::Internal& limitedAlphal
) const
{
    return pow
    (
        ((4*constant::mathematical::pi*n_)/3)
       *limitedAlphal/(1.0 + alphaNuc() - limitedAlphal),
        1.0/3.0
    );
}


Foam::dimensionedScalar
Foam::cavitationModels::SchnerrSauer::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::cavitationModels::SchnerrSauer::pCoeff
(
    const volScalarField::Internal& p
) const
{
    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal rho
    (
        limitedAlphal*rhol()
      + (scalar(1) - limitedAlphal)*rhov()
    );

    return
        (3*rhol()*rhov())*sqrt(2/(3*rhol()))
       *rRb(limitedAlphal)/(rho*sqrt(mag(p - pSat()) + 0.01*pSat()));
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::cavitationModels::SchnerrSauer::mDotcvAlpha() const
{
    const volScalarField::Internal& p =
        phases_.mesh().lookupObject<volScalarField>("p");

    const volScalarField::Internal pCoeff(this->pCoeff(p));

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField::Internal>>
    (
        Cc_*limitedAlphal*pCoeff*max(p - pSat(), p0_),
       -Cv_*(1.0 + alphaNuc() - limitedAlphal)*pCoeff*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::cavitationModels::SchnerrSauer::mDotcvP() const
{
    const volScalarField::Internal& p =
        phases_.mesh().lookupObject<volScalarField>("p");

    const volScalarField::Internal pCoeff(this->pCoeff(p));

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal apCoeff(limitedAlphal*pCoeff);

    return Pair<tmp<volScalarField::Internal>>
    (
        Cc_*(1.0 - limitedAlphal)*pos0(p - pSat())*apCoeff,
        (-Cv_)*(1.0 + alphaNuc() - limitedAlphal)*neg(p - pSat())*apCoeff
    );
}


void Foam::cavitationModels::SchnerrSauer::correct()
{}


bool Foam::cavitationModels::SchnerrSauer::read(const dictionary& dict)
{
    if (cavitationModel::read(dict))
    {
        dict.lookup("n") >> n_;
        dict.lookup("dNuc") >> dNuc_;
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

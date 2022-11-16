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
    const compressibleTwoPhases& phases
)
:
    cavitationModel(dict, phases),

    UInf_("UInf", dimVelocity, dict),
    tInf_("tInf", dimTime, dict),
    Cc_("Cc", dimless, dict),
    Cv_("Cv", dimless, dict),

    p0_("0", dimPressure, 0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::compressible::cavitationModels::Kunz::mDotcvAlphal() const
{
    const volScalarField::Internal& p =
        phases_.mesh().lookupObject<volScalarField>("p");

    const volScalarField::Internal mcCoeff_(Cc_*rhov()/tInf_);
    const volScalarField::Internal mvCoeff_
    (
        Cv_*rhov()/(0.5*rhol()*sqr(UInf_)*tInf_)
    );

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal pSatv(this->pSatv());
    const volScalarField::Internal pSatl(this->pSatl());

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff_*sqr(limitedAlphal)
       *max(p - pSatv, p0_)/max(p - pSatv, 0.01*pSatv),
       -mvCoeff_*min(p - pSatl, p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::compressible::cavitationModels::Kunz::mDotcvP() const
{
    const volScalarField::Internal& p =
        phases_.mesh().lookupObject<volScalarField>("p");

    const volScalarField::Internal mcCoeff_(Cc_*rhov()/tInf_);
    const volScalarField::Internal mvCoeff_
    (
        Cv_*rhov()/(0.5*rhol()*sqr(UInf_)*tInf_)
    );

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal pSatv(this->pSatv());
    const volScalarField::Internal pSatl(this->pSatl());

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff_*sqr(limitedAlphal)*(1 - limitedAlphal)
       *pos0(p - pSatv)/max(p - pSatv, 0.01*pSatv),
       -mvCoeff_*limitedAlphal*neg(p - pSatl)
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

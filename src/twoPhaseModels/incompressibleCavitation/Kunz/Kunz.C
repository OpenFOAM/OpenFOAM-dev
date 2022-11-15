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
namespace cavitationModels
{
    defineTypeNameAndDebug(Kunz, 0);
    addToRunTimeSelectionTable(cavitationModel, Kunz, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cavitationModels::Kunz::Kunz
(
    const dictionary& dict,
    const incompressibleTwoPhases& phases
)
:
    cavitationModel(dict, phases),

    UInf_("UInf", dimVelocity, dict),
    tInf_("tInf", dimTime, dict),
    Cc_("Cc", dimless, dict),
    Cv_("Cv", dimless, dict),

    p0_("0", pSat().dimensions(), 0.0),

    mcCoeff_(Cc_*rhov()/tInf_),
    mvCoeff_(Cv_*rhov()/(0.5*rhol()*sqr(UInf_)*tInf_))
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::cavitationModels::Kunz::mDotcvAlpha() const
{
    const volScalarField::Internal& p =
        phases_.mesh().lookupObject<volScalarField>("p");

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff_*sqr(limitedAlphal)
       *max(p - pSat(), p0_)/max(p - pSat(), 0.01*pSat()),
       -mvCoeff_*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::cavitationModels::Kunz::mDotcvP() const
{
    const volScalarField::Internal& p =
        phases_.mesh().lookupObject<volScalarField>("p");

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff_*sqr(limitedAlphal)*(1.0 - limitedAlphal)
       *pos0(p - pSat())/max(p - pSat(), 0.01*pSat()),
        (-mvCoeff_)*limitedAlphal*neg(p - pSat())
    );
}


void Foam::cavitationModels::Kunz::correct()
{}


bool Foam::cavitationModels::Kunz::read(const dictionary& dict)
{
    if (cavitationModel::read(dict))
    {
        dict.lookup("UInf") >> UInf_;
        dict.lookup("tInf") >> tInf_;
        dict.lookup("Cc") >> Cc_;
        dict.lookup("Cv") >> Cv_;

        mcCoeff_ = Cc_*rhov()/tInf_;
        mvCoeff_ = Cv_*rhov()/(0.5*rhol()*sqr(UInf_)*tInf_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

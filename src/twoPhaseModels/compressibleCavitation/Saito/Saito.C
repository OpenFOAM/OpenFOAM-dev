/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "Saito.H"
#include "addToRunTimeSelectionTable.H"

#include "mathematicalConstants.H"
using Foam::constant::mathematical::pi;

#include "physicoChemicalConstants.H"
using Foam::constant::physicoChemical::RR;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace cavitationModels
{
    defineTypeNameAndDebug(Saito, 0);
    addToRunTimeSelectionTable(cavitationModel, Saito, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressible::cavitationModels::Saito::Saito
(
    const dictionary& dict,
    const compressibleTwoPhases& phases
)
:
    cavitationModel(dict, phases),

    Ca_("Ca", dimless/dimLength, dict),
    Cv_("Cv", dimless, dict),
    Cc_("Cc", dimless, dict),
    alphaNuc_("alphaNuc", dimless, dict),

    p0_("0", dimPressure, 0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::cavitationModels::Saito::fT(const rhoThermo& thermo) const
{
    return sqrt(2*pi*(RR/thermo.W()()())*thermo.T()());
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::compressible::cavitationModels::Saito::mDotcvAlphal() const
{
    const volScalarField::Internal& p = thermol().p();

    const volScalarField::Internal alphav
    (
        min(max(this->alphav(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal alphal
    (
        min(max(this->alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal alphavNuc(max(alphav, alphaNuc_));

    const volScalarField::Internal A(Ca_*alphal*alphavNuc);

    const volScalarField::Internal mvCoeff(Cv_*A*rhol()/(rhov()*fT(thermol())));
    const volScalarField::Internal mcCoeff(Cc_*A/fT(thermol()));

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff*alphal*max(p - pSatv(), p0_),
       -mvCoeff*alphavNuc*min(p - pSatl(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::compressible::cavitationModels::Saito::mDotcvP() const
{
    const volScalarField::Internal& p = thermol().p();

    const volScalarField::Internal alphav
    (
        min(max(this->alphav(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal alphal
    (
        min(max(this->alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal alphavNuc(max(alphav, alphaNuc_));

    const volScalarField::Internal A(Ca_*alphal*alphavNuc);

    const volScalarField::Internal mvCoeff(Cv_*A*rhol()/(rhov()*fT(thermol())));
    const volScalarField::Internal mcCoeff(Cc_*A/fT(thermol()));

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff*alphal*alphav*pos0(p - pSatv()),
       -mvCoeff*alphal*alphavNuc*neg(p - pSatl())
    );
}


void Foam::compressible::cavitationModels::Saito::correct()
{}


bool Foam::compressible::cavitationModels::Saito::read
(
    const dictionary& dict
)
{
    if (cavitationModel::read(dict))
    {
        dict.lookup("Ca") >> Ca_;
        dict.lookup("Cv") >> Cv_;
        dict.lookup("Cc") >> Cc_;
        dict.lookup("alphaNuc") >> alphaNuc_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

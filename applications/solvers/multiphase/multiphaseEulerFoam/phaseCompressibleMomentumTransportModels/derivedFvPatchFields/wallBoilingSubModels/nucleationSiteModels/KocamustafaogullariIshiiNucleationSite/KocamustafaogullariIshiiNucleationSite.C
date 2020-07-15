/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "KocamustafaogullariIshiiNucleationSite.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleationSiteModels
{
    defineTypeNameAndDebug(KocamustafaogullariIshiiNucleationSite, 0);
    addToRunTimeSelectionTable
    (
        nucleationSiteModel,
        KocamustafaogullariIshiiNucleationSite,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::KocamustafaogullariIshiiNucleationSite
(
    const dictionary& dict
)
:
    nucleationSiteModel(),
    Cn_(dict.lookupOrDefault<scalar>("Cn", 1))
{}


Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::KocamustafaogullariIshiiNucleationSite
(
    const KocamustafaogullariIshiiNucleationSite& model
)
:
    nucleationSiteModel(),
    Cn_(model.Cn_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::
~KocamustafaogullariIshiiNucleationSite()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::N
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L,
    const scalarField& dDep,
    const scalarField& fDep
) const
{
    const fvPatchScalarField& Tw =
        liquid.thermo().T().boundaryField()[patchi];

    const scalarField rhoLiquid(liquid.thermo().rho(patchi));
    const scalarField rhoVapor(vapor.thermo().rho(patchi));
    const scalarField rhoM((rhoLiquid - rhoVapor)/rhoVapor);

    const scalarField sigmaw
    (
        liquid.fluid().sigma(phasePairKey(liquid.name(), vapor.name()), patchi)
    );

    //eq. (32)
    const scalarField f(2.157e-7*pow(rhoM,-3.2)*pow(1 + 0.0049*rhoM,4.13));

    // eq. (17)
    const scalarField rRc(max(Tw-Tsatw,scalar(0))*rhoVapor*L/(2*sigmaw*Tsatw));

    return (Cn_/sqr(dDep))*pow(rRc,4.4)*f;
}


void Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::write(Ostream& os) const
{
    nucleationSiteModel::write(os);
    writeKeyword(os, "Cn") << Cn_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //

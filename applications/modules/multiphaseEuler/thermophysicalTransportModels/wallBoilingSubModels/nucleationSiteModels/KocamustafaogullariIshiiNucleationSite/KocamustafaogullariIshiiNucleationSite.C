/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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
#include "wallBoilingModelsCoefficient.H"
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ScalarFieldType>
Foam::tmp<ScalarFieldType>
Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::calculate
(
    const ScalarFieldType& Tsatw,
    const ScalarFieldType& L,
    const ScalarFieldType& dDep,
    const ScalarFieldType& Tw,
    const ScalarFieldType& rhoLiquid,
    const ScalarFieldType& rhoVapour,
    const ScalarFieldType& sigma
) const
{
    const ScalarFieldType rhoM((rhoLiquid - rhoVapour)/rhoVapour);

    const dimensionedScalar zeroT_(dimTemperature, 0);
    auto zeroT = coefficient<ScalarFieldType>::value(zeroT_);

    auto Cn = coefficient<ScalarFieldType>::value(Cn_);

    // eq. (32)
    const ScalarFieldType f
    (
        2.157e-7*pow(rhoM, -3.2)*pow(1 + 0.0049*rhoM, 4.13)
    );

    // eq. (17)
    const ScalarFieldType rRc
    (
        dDep*max(Tw - Tsatw, zeroT)*rhoVapour*L/(4*sigma*Tsatw)
    );

    return Cn/sqr(dDep)*pow(rRc, 4.4)*f;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::KocamustafaogullariIshiiNucleationSite
(
    const dictionary& dict
)
:
    nucleationSiteModel(),
    Cn_("Cn", dimless, dict, 1)
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
KocamustafaogullariIshiiNucleationSite::nucleationSiteDensity
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L,
    const scalarField& dDep,
    const scalarField& fDep
) const
{
    const scalarField& Tw =
        liquid.thermo().T().boundaryField()[patchi];

    return
        calculate
        (
            Tsatw,
            L,
            dDep,
            Tw,
            static_cast<const scalarField&>
            (
                liquid.rho().boundaryField()[patchi]
            ),
            static_cast<const scalarField&>
            (
                vapour.rho().boundaryField()[patchi]
            ),
            liquid.fluid().sigma(phaseInterfaceKey(liquid, vapour), patchi)()
        );
}


Foam::tmp<Foam::volScalarField>
Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::nucleationSiteDensity
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const phaseModel& solid,
    const volScalarField& Tf,
    const volScalarField& Tsatw,
    const volScalarField& L,
    const volScalarField& dDep,
    const volScalarField& fDep
) const
{
    return calculate
    (
        Tsatw,
        L,
        dDep,
        Tf,
        liquid.rho(),
        vapour.rho(),
        liquid.fluid().sigma(phaseInterfaceKey(liquid, vapour))()
    );
}


void Foam::wallBoilingModels::nucleationSiteModels::
KocamustafaogullariIshiiNucleationSite::write(Ostream& os) const
{
    nucleationSiteModel::write(os);
    writeKeyword(os, "Cn") << Cn_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //

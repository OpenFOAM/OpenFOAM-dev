/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "LemmertChawla.H"
#include "addToRunTimeSelectionTable.H"
#include "wallBoilingModelsCoefficient.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleationSiteModels
{
    defineTypeNameAndDebug(LemmertChawla, 0);
    addToRunTimeSelectionTable
    (
        nucleationSiteModel,
        LemmertChawla,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ScalarFieldType>
Foam::tmp<ScalarFieldType>
Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::calculate
(
    const ScalarFieldType& Tsatw,
    const ScalarFieldType& Tw
) const
{
    auto Cn = coefficient<ScalarFieldType>::value(Cn_);
    auto NRef = coefficient<ScalarFieldType>::value(NRef_);
    auto deltaTRef = coefficient<ScalarFieldType>::value(deltaTRef_);

    return Cn*NRef*pow(max((Tw - Tsatw)/deltaTRef, scalar(0)), 1.805);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::LemmertChawla
(
    const dictionary& dict
)
:
    nucleationSiteModel(),
    Cn_(dimensionedScalar::lookupOrDefault("Cn", dict, dimless, 1)),
    NRef_
    (
        dimensionedScalar::lookupOrDefault
        (
            "NRef",
            dict,
            dimless/dimArea,
            9.922e5
        )
    ),
    deltaTRef_
    (
        dimensionedScalar::lookupOrDefault
        (
            "deltaTRef",
            dict,
            dimTemperature,
            10
        )
    )
{}


Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::LemmertChawla
(
    const LemmertChawla& model
)
:
    nucleationSiteModel(),
    Cn_(model.Cn_),
    NRef_(model.NRef_),
    deltaTRef_(model.deltaTRef_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::~LemmertChawla()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::
nucleationSiteDensity
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

    return calculate(Tsatw, Tw);
}


Foam::tmp<Foam::volScalarField>
Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::
nucleationSiteDensity
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
    return calculate(Tsatw, Tf);
}


void Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::write
(
    Ostream& os
) const
{
    nucleationSiteModel::write(os);
    writeKeyword(os, "Cn") << Cn_ << token::END_STATEMENT << nl;
    writeKeyword(os, "NRef") << NRef_ << token::END_STATEMENT << nl;
    writeKeyword(os, "deltaTRef") << deltaTRef_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //

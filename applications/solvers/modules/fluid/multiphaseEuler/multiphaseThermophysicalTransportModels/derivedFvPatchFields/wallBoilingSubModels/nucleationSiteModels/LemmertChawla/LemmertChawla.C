/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::LemmertChawla
(
    const dictionary& dict
)
:
    nucleationSiteModel(),
    Cn_(dict.lookupOrDefault<scalar>("Cn", 1)),
    NRef_(dict.lookupOrDefault<scalar>("NRef", 9.922e5)),
    deltaTRef_(dict.lookupOrDefault<scalar>("deltaTRef", 10))
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
Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::N
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

    return Cn_*NRef_*pow(max((Tw - Tsatw)/deltaTRef_, scalar(0)), 1.805);
}


void Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla::
    write(Ostream& os) const
{
    nucleationSiteModel::write(os);
    writeKeyword(os, "Cn") << Cn_ << token::END_STATEMENT << nl;
    writeKeyword(os, "NRef") << NRef_ << token::END_STATEMENT << nl;
    writeKeyword(os, "deltaTRef") << deltaTRef_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //

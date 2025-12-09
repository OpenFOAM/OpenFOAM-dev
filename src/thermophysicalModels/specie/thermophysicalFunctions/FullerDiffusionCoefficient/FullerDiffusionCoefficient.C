/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "FullerDiffusionCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{
    addScalarFunction2(FullerDiffusionCoefficient);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function2s::FullerDiffusionCoefficient::FullerDiffusionCoefficient
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction2<scalar, FullerDiffusionCoefficient>(name),
    Va_(dict.lookup<scalar>("Va")),
    Vb_(dict.lookup<scalar>("Vb")),
    Ma_(dict.lookup<scalar>("Ma")),
    Mb_(dict.lookup<scalar>("Mb")),
    alpha_(sqrt(1/Ma_ + 1/Mb_)),
    beta_(sqr((cbrt(Va_) + cbrt(Vb_))))
{
    assertNoConvertUnits(typeName, units, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function2s::FullerDiffusionCoefficient::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, "Va", Va_);
    writeEntry(os, "Vb", Vb_);
    writeEntry(os, "Ma", Ma_);
    writeEntry(os, "Mb", Mb_);
}


// ************************************************************************* //

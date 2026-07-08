/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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

#include "Antoine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    addScalarFunction1(Antoine);
    addScalarFunction1(logAntoine);
    addScalarFunction1(inverseAntoine);
}
}


const Foam::scalar Foam::Function1s::AntoineCoeffs::ln10_ = log(scalar(10));


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::AntoineCoeffs::AntoineCoeffs
(
    const word& typeName,
    const unitSets& units,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A", dimless)),
    B_(dict.lookup<scalar>("B", units.x)),
    C_(dict.lookup<scalar>("C", units.x))
{
    assertNoConvertUnits(typeName, {units::any, units.value}, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function1s::AntoineCoeffs::write
(
    Ostream& os,
    const unitSets& units
) const
{
    writeEntry(os, "A", units::unitless, A_);
    writeEntry(os, "B", units.x, B_);
    writeEntry(os, "C", units.x, C_);
}


// ************************************************************************* //

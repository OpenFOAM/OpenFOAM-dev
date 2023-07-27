/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "reverseRamp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    addScalarFunction1(reverseRamp);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::reverseRamp::reverseRamp
(
    const word& name,
    const dictionary& dict
)
:
    Ramp<reverseRamp>(name, dict),
    ramp_(Function1<scalar>::New("ramp", dict))
{}


Foam::Function1s::reverseRamp::reverseRamp
(
    const reverseRamp& rr
)
:
    Ramp<reverseRamp>(rr),
    ramp_(rr.ramp_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1s::reverseRamp::~reverseRamp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function1s::reverseRamp::write(Ostream& os) const
{
    Ramp<reverseRamp>::write(os);
    writeEntry(os, ramp_());
}


// ************************************************************************* //

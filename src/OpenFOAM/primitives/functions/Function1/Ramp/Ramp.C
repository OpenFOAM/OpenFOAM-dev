/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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

#include "Ramp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Function1Type>
Foam::Function1s::Ramp<Function1Type>::Ramp
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<scalar, Function1Type>(name),
    start_(dict.lookupOrDefault<scalar>("start", units.x, 0)),
    duration_(dict.lookup<scalar>("duration", units.x))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Function1Type>
Foam::Function1s::Ramp<Function1Type>::~Ramp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Function1Type>
void Foam::Function1s::Ramp<Function1Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, "start", units.x, start_);
    writeEntry(os, "duration", units.x, duration_);
}


// ************************************************************************* //

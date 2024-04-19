/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "Square.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Square<Type>::Square
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Square<Type>>(name),
    amplitude_(Function1<Type>::New("amplitude", units, dict)),
    frequency_(dict.lookup<scalar>("frequency", unitless/units.x)),
    start_(dict.lookupOrDefault<scalar>("start", units.x, 0)),
    level_(Function1<Type>::New("level", units, dict)),
    markSpace_(dict.lookupOrDefault<scalar>("markSpace", unitless, 1)),
    integrable_(amplitude_->constant() && level_->constant())
{}


template<class Type>
Foam::Function1s::Square<Type>::Square(const Square<Type>& se)
:
    FieldFunction1<Type, Square<Type>>(se),
    amplitude_(se.amplitude_, false),
    frequency_(se.frequency_),
    start_(se.start_),
    level_(se.level_, false),
    markSpace_(se.markSpace_),
    integrable_(se.integrable_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Square<Type>::~Square()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Square<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, units, amplitude_());
    writeEntry(os, "frequency", unitless/units.x, frequency_);
    writeEntry(os, "start", units.x, start_);
    writeEntry(os, units, level_());
    writeEntry(os, "markSpace", unitless, markSpace_);
}


// ************************************************************************* //

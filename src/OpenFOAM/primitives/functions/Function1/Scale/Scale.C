/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "Scale.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Scale<Type>::Scale
(
    const word& name,
    const Function1<scalar>& scale,
    const Function1<scalar>& xScale,
    const Function1<Type>& value
)
:
    FieldFunction1<Type, Scale<Type>>(name),
    scale_(scale.clone().ptr()),
    constantScale_(scale_->constant()),
    xScale_(xScale.clone().ptr()),
    constantXScale_(xScale_->constant()),
    value_(value.clone().ptr()),
    constantValue_(value_->constant())
{}


template<class Type>
Foam::Function1s::Scale<Type>::Scale
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Scale<Type>>(name),
    scale_(Function1<scalar>::New("scale", units.x, unitAny, dict)),
    constantScale_(scale_->constant()),
    xScale_
    (
        dict.found("xScale")
      ? Function1<scalar>::New("xScale", units.x, unitless, dict)
      : autoPtr<Function1<scalar>>(new Constant<scalar>("xScale", 1))
    ),
    constantXScale_(xScale_->constant()),
    value_(Function1<Type>::New("value", units.x, unitAny, dict)),
    constantValue_(value_->constant())
{}


template<class Type>
Foam::Function1s::Scale<Type>::Scale(const Scale<Type>& se)
:
    FieldFunction1<Type, Scale<Type>>(se),
    scale_(se.scale_, false),
    constantScale_(se.constantScale_),
    xScale_(se.xScale_, false),
    constantXScale_(se.constantXScale_),
    value_(se.value_, false),
    constantValue_(se.constantValue_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Scale<Type>::~Scale()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Scale<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, units.x, unitless, scale_());
    writeEntry(os, units.x, units.x, xScale_());
    writeEntry(os, units, value_());
}


// ************************************************************************* //

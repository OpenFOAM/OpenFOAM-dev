/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "Scale2.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Function2s::Scale<Type>::read(const dictionary& dict)
{
    scale_ = Function2<scalar>::New("scale", dict);
    xScale_ =
        dict.found("xScale")
      ? Function1<scalar>::New("xScale", dict)
      : autoPtr<Function1<scalar>>
        (
            new Function1s::Constant<scalar>("xScale", 1)
        );
    yScale_ =
        dict.found("yScale")
      ? Function1<scalar>::New("yScale", dict)
      : autoPtr<Function1<scalar>>
        (
            new Function1s::Constant<scalar>("yScale", 1)
        );
    value_ = Function2<Type>::New("value", dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Scale<Type>::Scale
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction2<Type, Scale<Type>>(name)
{
    read(dict);
}


template<class Type>
Foam::Function2s::Scale<Type>::Scale(const Scale<Type>& se)
:
    FieldFunction2<Type, Scale<Type>>(se),
    scale_(se.scale_, false),
    xScale_(se.xScale_, false),
    yScale_(se.yScale_, false),
    value_(se.value_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Scale<Type>::~Scale()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function2s::Scale<Type>::write(Ostream& os) const
{
    writeEntry(os, scale_());
    writeEntry(os, xScale_());
    writeEntry(os, yScale_());
    writeEntry(os, value_());
}


// ************************************************************************* //

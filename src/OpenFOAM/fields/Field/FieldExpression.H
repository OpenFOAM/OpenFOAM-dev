/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "expressions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define Field to be expressionable by iteration
template<class Type>
struct Expressionable<Field<Type>>
:
    public ExpressionableIterable<Field<Type>>
{};

//- Define Field operations to be expressionable by iteration
template<class Type, class ... Args>
struct Expressionable<Field<Type>, Arguments<Args ...>>
:
    public ExpressionableIterable<Field<Type>, Arguments<Args ...>>
{};

//- Define Field to operate by iteration
template<class Type, class ... Args>
struct Operate<Field<Type>, Arguments<Args ...>>
:
    public OperateIterable<Field<Type>, Arguments<Args ...>>
{};

//- Define standard iterable types returned by Field expressions
template<class Type>
struct Result<Field<Type>>
:
    public ResultIterable<Field<Type>>
{};

//- Define iterable reduction of a Field
template<class Type>
struct Reduce<Field<Type>>
:
    public ReduceIterable<Field<Type>>
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //

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

//- Define PtrField to be expressionable by iteration
template<class Type>
struct Expressionable<PtrField<Type>>
:
    public ExpressionableIterable<PtrField<Type>>
{};

//- Define PtrField operations to be expressionable by iteration
template<class Type, class ... Args>
struct Expressionable<PtrField<Type>, Arguments<Args ...>>
:
    public ExpressionableIterable<PtrField<Type>, Arguments<Args ...>>
{};

//- Define PtrField to operate by iteration
template<class Type, class ... Args>
struct Operate<PtrField<Type>, Arguments<Args ...>>
:
    public OperateIterable<PtrField<Type>, Arguments<Args ...>>
{};

//- Define standard iterable types returned by PtrField expressions
template<class Type>
struct Result<PtrField<Type>>
:
    public ResultIterable<PtrField<Type>>
{};

//- Define iterable reduction of a PtrField
template<class Type>
struct Reduce<PtrField<Type>>
:
    public ReduceIterable<PtrField<Type>>
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //

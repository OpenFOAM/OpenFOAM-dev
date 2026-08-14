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

//- Define SubField to be expressionable by iteration
template<class Type>
struct Expressionable<SubField<Type>>
:
    public ExpressionableIterable<SubField<Type>>
{};

//- Define SubField operations to be expressionable by iteration
template<class Type, class ... Args>
struct Expressionable<SubField<Type>, Arguments<Args ...>>
:
    public ExpressionableIterable<SubField<Type>, Arguments<Args ...>>
{};

//- Define SubField to operate by iteration
template<class Type, class ... Args>
struct Operate<SubField<Type>, Arguments<Args ...>>
:
    public OperateIterable<SubField<Type>, Arguments<Args ...>>
{};

//- Define the types returned by SubField expressions as the same as for
//  Field expressions
template<class Type>
struct Result<SubField<Type>>
:
    public Result<Field<Type>>
{};

//- Define reduction of SubField expressions as the same as for Field
//  expressions
template<class Type>
struct Reduce<SubField<Type>>
:
    public Reduce<Field<Type>>
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //

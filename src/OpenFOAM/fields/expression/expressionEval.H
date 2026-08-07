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

Namespace
    Foam

Description
    A function which evaluates an expression into a temporary field. Note that
    this is in the top-level namespace, not the expression sub-namespace. This
    is because this is designed to be used in high-level code, so it is
    desirable to avoid the clutter associated with the sub-namespace.

\*---------------------------------------------------------------------------*/

#ifndef expressionEval_H
#define expressionEval_H

#include "expressionContainers.H"
#include "expressionResult.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Evaluate an expression and return a temporary object
template<class Expression, class = EnableIfExpressionable<Expression>>
tmp
<
    typename expression::Result
    <
        typename expression::AllContainers<Expression>::type
    >::template type<Expression>
>
eval(const Expression& e)
{
    return
        tmp
        <
            typename expression::Result
            <
                typename expression::AllContainers<Expression>::type
            >::template type<Expression>
        >
        (
            new typename expression::Result
            <
                typename expression::AllContainers<Expression>::type
            >::template type<Expression>(e)
        );
}


//- Evaluate an expression and write to Ostream
template<class Expression, class = EnableIfExpressionable<Expression>>
Ostream& operator<<(Ostream& os, const Expression& e)
{
    os << eval(e);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

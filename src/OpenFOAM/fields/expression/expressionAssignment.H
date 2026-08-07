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

Class
    Foam::expression::Assignment

Description
    Simple wrapper that holds an reference to an expression and permits an
    expressionable container to assign to it. Adds virtualisation so that the
    container's assignment operator may be virtual.

\*---------------------------------------------------------------------------*/

#ifndef expressionAssignment_H
#define expressionAssignment_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                          Class Assignment Declaration
\*---------------------------------------------------------------------------*/

template<class Container>
class Assignment
{
public:

    virtual ~Assignment()
    {}

    virtual void assign(Container& c) const = 0;
};

/*---------------------------------------------------------------------------*\
                      Class AssignmentToExpression Declaration
\*---------------------------------------------------------------------------*/

template<class Container, class Expression>
class AssignmentToExpression
:
    public Assignment<Container>
{
    const Expression& e_;

public:

    AssignmentToExpression(const Expression& e)
    :
        e_(e)
    {}

    virtual void assign(Container& c) const
    {
        c = e_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return an assignment object for an expression
template<class Container, class Expression>
AssignmentToExpression<Container, Expression> assignment(const Expression& e)
{
    return AssignmentToExpression<Container, Expression>(e);
}

//- Assign a container to an assignment. Implemented as the '==' operator as
//  C++ does not permit non-member overloads of '='.
template<class Container>
void operator==(Container& c, const Assignment<Container>& a)
{
    a.assign(c);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

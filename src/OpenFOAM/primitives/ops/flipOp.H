/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::flipOp

Description
    Class containing functor to negate primitives. Dummy for all other types.

    Used in mesh transformations where face can flip.

SourceFiles
    flipOp.C

\*---------------------------------------------------------------------------*/

#ifndef flipOp_H
#define flipOp_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flipOp Declaration
\*---------------------------------------------------------------------------*/

class flipOp
{
public:

    template<class Type>
    Type operator()(const Type& val) const
    {
        return val;
    }
};


class noOp
{
public:

    template<class Type>
    Type operator()(const Type& val) const
    {
        return val;
    }
};


class flipLabelOp
{
public:

    label operator()(const label& val) const
    {
        return -val;
    }
};


// Template specialisation for primitives that support negation
template<> scalar flipOp::operator()(const scalar&) const;
template<> vector flipOp::operator()(const vector&) const;
template<> sphericalTensor flipOp::operator()(const sphericalTensor&) const;
template<> symmTensor flipOp::operator()(const symmTensor&) const;
template<> tensor flipOp::operator()(const tensor&) const;
template<> triad flipOp::operator()(const triad&) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

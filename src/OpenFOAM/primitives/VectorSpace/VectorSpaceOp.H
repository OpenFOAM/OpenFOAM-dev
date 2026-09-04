/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::VectorSpaceOps

Description
    Operator function for VectorSpace

\*---------------------------------------------------------------------------*/

#ifndef VectorSpaceOp_H
#define VectorSpaceOp_H

#include "direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class VectorSpaceOp Declaration
\*---------------------------------------------------------------------------*/

template<direction N, direction I>
struct VectorSpaceOp
{
    template<class Op, class ... Vs>
    static inline void op(const Op& op, Vs& ... vs)
    {
        op(vs.v_[I] ...);
        VectorSpaceOp<N, I + 1>::op(op, vs ...);
    }
};


template<direction N>
struct VectorSpaceOp<N, N>
{
    template<class Op, class ... Vs>
    static inline void op(const Op& op, Vs& ... vs)
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

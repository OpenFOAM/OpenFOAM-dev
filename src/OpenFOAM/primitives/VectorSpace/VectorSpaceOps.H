/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Operator functions for VectorSpace.

\*---------------------------------------------------------------------------*/

#ifndef VectorSpaceOps_H
#define VectorSpaceOps_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<direction N, direction I>
class VectorSpaceOps
{
public:

    template<class V, class S, class EqOp>
    static inline void eqOpS(V& vs, const S& s, EqOp eo)
    {
        eo(vs.v_[I], s);
        VectorSpaceOps<N, I + 1>::eqOpS(vs, s, eo);
    }

    template<class S, class V, class EqOp>
    static inline void SeqOp(S& s, const V& vs, EqOp eo)
    {
        eo(s, vs.v_[I]);
        VectorSpaceOps<N, I + 1>::SeqOp(s, vs, eo);
    }

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1& vs1, const V2& vs2, EqOp eo)
    {
        eo(vs1.v_[I], vs2.v_[I]);
        VectorSpaceOps<N, I + 1>::eqOp(vs1, vs2, eo);
    }


    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1& vs1, const S& s, Op o)
    {
        vs.v_[I] = o(vs1.v_[I], s);
        VectorSpaceOps<N, I + 1>::opVS(vs, vs1, s, o);
    }

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S& s,  const V1& vs1, Op o)
    {
        vs.v_[I] = o(s, vs1.v_[I]);
        VectorSpaceOps<N, I + 1>::opSV(vs, s, vs1, o);
    }

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1& vs1, const V1& vs2, Op o)
    {
        vs.v_[I] = o(vs1.v_[I], vs2.v_[I]);
        VectorSpaceOps<N, I + 1>::op(vs, vs1, vs2, o);
    }
};


template<direction N>
class VectorSpaceOps<N, N>
{
public:

    template<class V, class S, class EqOp>
    static inline void eqOpS(V&, const S&, EqOp)
    {}

    template<class S, class V, class EqOp>
    static inline void SeqOp(S&, const V&, EqOp)
    {}

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1&, const V2&, EqOp)
    {}


    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1&, const S&, Op)
    {}

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S&, const V1&, Op)
    {}

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1&, const V1&, Op)
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

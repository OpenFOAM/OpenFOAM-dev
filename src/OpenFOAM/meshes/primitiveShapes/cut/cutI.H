/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

#include "cut.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Modify a uniform operation for reordering a tri (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& triReorder
(
    const cut::uniformOp<Type>& x,
    const FixedList<label, 3>&
)
{
    return x;
}


//- Modify a uniform operation for cutting a tri from a tri (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& triCutTri
(
    const cut::uniformOp<Type>& x,
    const Pair<scalar>&
)
{
    return x;
}


//- Modify a uniform operation for cutting a quad from a tri (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& triCutQuad
(
    const cut::uniformOp<Type>& x,
    const Pair<scalar>&
)
{
    return x;
}


//- Modify a uniform operation for reordering a tet (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& tetReorder
(
    const cut::uniformOp<Type>& x,
    const FixedList<label, 4>&
)
{
    return x;
}


//- Modify a uniform operation for cutting a tet from a tet (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& tetCutTet
(
    const cut::uniformOp<Type>& x,
    const FixedList<scalar, 3>&
)
{
    return x;
}


//- Modify a uniform operation for cutting prism0 from a tet (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& tetCutPrism0
(
    const cut::uniformOp<Type>& x,
    const FixedList<scalar, 3>&
)
{
    return x;
}


//- Modify a uniform operation for cutting prism01 from a tet (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& tetCutPrism01
(
    const cut::uniformOp<Type>& x,
    const FixedList<scalar, 4>&
)
{
    return x;
}


//- Modify a uniform operation for cutting prism23 from a tet (does nothing)
template<class Type>
inline const cut::uniformOp<Type>& tetCutPrism23
(
    const cut::uniformOp<Type>& x,
    const FixedList<scalar, 4>&
)
{
    return x;
}


//- Modify a fixed list for reordering a tri (does nothing)
template<class Type, unsigned Size>
inline FixedList<Type, 3> triReorder
(
    const FixedList<Type, Size>& x,
    const FixedList<label, 3>& indices
)
{
    FixedList<Type, 3> result;
    for (unsigned i = 0; i < 3; ++ i)
    {
        result[i] = x[indices[i]];
    }
    return result;
}


//- Modify a list for cutting a tri from a tri
template<class Type>
inline FixedList<Type, 3> triCutTri
(
    const FixedList<Type, 3>& x,
    const Pair<scalar>& f
)
{
    FixedList<Type, 3> result;
    result[0] = x[0];
    for (label i = 0; i < 2; ++ i)
    {
        result[i+1] = (1 - f[i])*x[0] + f[i]*x[i+1];
    }
    return result;
}


//- Modify a list for cutting a quad from a tri
template<class Type>
inline FixedList<Type, 4> triCutQuad
(
    const FixedList<Type, 3>& x,
    const Pair<scalar>& f
)
{
    FixedList<Type, 4> result;
    for (label i = 0; i < 2; ++ i)
    {
        result[i] = x[i+1];
        result[3-i] = (1 - f[i])*x[0] + f[i]*x[i+1];
    }
    return result;
}


//- Modify a fixed list for reordering a tet (does nothing)
template<class Type, unsigned Size>
inline FixedList<Type, 4> tetReorder
(
    const FixedList<Type, Size>& x,
    const FixedList<label, 4>& indices
)
{
    FixedList<Type, 4> result;
    for (unsigned i = 0; i < 4; ++ i)
    {
        result[i] = x[indices[i]];
    }
    return result;
}


//- Modify a list for cutting a tet from a tet
template<class Type>
inline FixedList<Type, 4> tetCutTet
(
    const FixedList<Type, 4>& x,
    const FixedList<scalar, 3>& f
)
{
    FixedList<Type, 4> result;
    result[0] = x[0];
    for (label i = 0; i < 3; ++ i)
    {
        result[i+1] = (1 - f[i])*x[0] + f[i]*x[i+1];
    }
    return result;
}


//- Modify a list for cutting prism0 from a tet
template<class Type>
inline FixedList<Type, 6> tetCutPrism0
(
    const FixedList<Type, 4>& x,
    const FixedList<scalar, 3>& f
)
{
    FixedList<Type, 6> result;
    for (label i = 0; i < 3; ++ i)
    {
        result[i] = (1 - f[i])*x[0] + f[i]*x[i+1];
        result[i+3] = x[i+1];
    }
    return result;
}


//- Modify a list for cutting prism01 from a tet
template<class Type>
inline FixedList<Type, 6> tetCutPrism01
(
    const FixedList<Type, 4>& x,
    const FixedList<scalar, 4>& f
)
{
    FixedList<Type, 6> result;
    for (label i = 0; i < 2; ++ i)
    {
        result[3*i] = x[i];
        for (label j = 0; j < 2; ++ j)
        {
            result[3*i+j+1] = (1 - f[2*i+j])*x[i] + f[2*i+j]*x[j+2];
        }
    }

    return result;
}


//- Modify a list for cutting prism23 from a tet
template<class Type>
inline FixedList<Type, 6> tetCutPrism23
(
    const FixedList<Type, 4>& x,
    const FixedList<scalar, 4>& f
)
{
    FixedList<Type, 6> result = tetCutPrism01(x, f);
    result[0] = x[2];
    result[3] = x[3];
    Swap(result[2], result[4]);
    return result;
}


//- Cut a tri from a tri and apply an operation to the result. The cut is made
//  along the two edges connected to vertex 0, and the cut locations are given
//  as factors along these edges. The result is the side connected to vertex 0.
template<class Op, class Point>
inline typename Op::result triCutTri
(
    const Op& op,
    const FixedList<Point, 3>& p,
    const Pair<scalar>& f
)
{
    return Op(triCutTri(op, f))(triCutTri(p, f));
}


//- Apply an operation to a quad. Splits the quad into two tris.
template<class Op, class OpData, class Point>
inline typename Op::result quadOp
(
    const OpData& opData,
    const FixedList<Point, 4>& p
)
{
    static const FixedList<FixedList<label, 3>, 2> i =
        {{0, 1, 2}, {0, 2, 3}};
    return
        Op(triReorder(opData, i[0]))(triReorder(p, i[0]))
      + Op(triReorder(opData, i[1]))(triReorder(p, i[1]));
}


//- Cut a quad from a tri and apply an operation to the result. The cuts are
//  the same as for triCutTri. The result is the side connected to vertices 1
//  and 2.
template<class Op, class Point>
inline typename Op::result triCutQuad
(
    const Op& op,
    const FixedList<Point, 3>& p,
    const FixedList<scalar, 2>& f
)
{
    return quadOp<Op>(triCutQuad(op, f), triCutQuad(p, f));
}


//- Cut a tet from a tet and apply an operation to the result. The cut is made
//  along the three edges connected to vertex 0, and the cut locations are given
//  as factors along these edges. The result is the side connected to vertex 0.
template<class Op, class Point>
inline typename Op::result tetCutTet
(
    const Op& op,
    const FixedList<Point, 4>& p,
    const FixedList<scalar, 3>& f
)
{
    return Op(tetCutTet(op, f))(tetCutTet(p, f));
}


//- Apply an operation to a prism. Splits the prism into three tets.
template<class Op, class OpData, class Point>
inline typename Op::result prismOp
(
    const OpData& opData,
    const FixedList<Point, 6>& p
)
{
    static const FixedList<FixedList<label, 4>, 3> i =
        {{0, 1, 2, 4}, {0, 2, 5, 4}, {0, 4, 5, 3}};
    return
        Op(tetReorder(opData, i[0]))(tetReorder(p, i[0]))
      + Op(tetReorder(opData, i[1]))(tetReorder(p, i[1]))
      + Op(tetReorder(opData, i[2]))(tetReorder(p, i[2]));
}


//- Cut a prism from a tet and apply an operation to the result. The cuts are
//  the same as for tetCutTet. The result is the side connected to vertices 1,
//  2 and 3.
template<class Op, class Point>
inline typename Op::result tetCutPrism0
(
    const Op& op,
    const FixedList<Point, 4>& p,
    const FixedList<scalar, 3>& f
)
{
    return prismOp<Op>(tetCutPrism0(op, f), tetCutPrism0(p, f));
}


//- Cut a prism from a tet and apply an operation to the result. The cut is made
//  along four edges, not edges 01 or 23, and the cut locations are given as
//  factors along these edges. The result is the side connected to edge 01.
template<class Op, class Point>
inline typename Op::result tetCutPrism01
(
    const Op& op,
    const FixedList<Point, 4>& p,
    const FixedList<scalar, 4>& f
)
{
    return prismOp<Op>(tetCutPrism01(op, f), tetCutPrism01(p, f));
}


//- Cut a prism from a tet and apply an operation to the result. The cuts are
//  the same as for tetCutPrism01. The result is the side connected to edge 23.
template<class Op, class Point>
inline typename Op::result tetCutPrism23
(
    const Op& op,
    const FixedList<Point, 4>& p,
    const FixedList<scalar, 4>& f
)
{
    return prismOp<Op>(tetCutPrism23(op, f), tetCutPrism23(p, f));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

template<class Point, class AboveOp, class BelowOp>
typename Foam::cut::opAddResult<AboveOp, BelowOp>::type Foam::triCut
(
    const FixedList<Point, 3>& tri,
    const FixedList<scalar, 3>& level,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
)
{
    // If everything is positive or negative, then process the triangle as a
    // whole, and do a quick return
    if (level[0] >= 0 && level[1] >= 0 && level[2] >= 0)
    {
        return aboveOp(tri) + belowOp();
    }
    if (level[0] <= 0 && level[1] <= 0 && level[2] <= 0)
    {
        return aboveOp() + belowOp(tri);
    }

    // There will be just one edge without a sign change. Find it, and put it
    // opposite the first vertex. This may change the sign of the tri.
    FixedList<label, 3> indices({0, 1, 2});
    label i;
    for (i = 0; i < 3; ++ i)
    {
        if (level[(i + 1)%3]*level[(i + 2)%3] >= 0)
        {
            Swap(indices[0], indices[i]);
            break;
        }
    }
    if (i == 3)
    {
        FatalErrorInFunction
            << "The number of tri vertices above the level set should always "
            << "be 1" << exit(FatalError);
    }

    // Correct the sign
    if (indices[0] != 0)
    {
        Swap(indices[1], indices[2]);
    }

    // Permute the data
    const FixedList<Point, 3> p = triReorder(tri, indices);
    const FixedList<scalar, 3> l = triReorder(level, indices);
    AboveOp a = triReorder(aboveOp, indices);
    BelowOp b = triReorder(belowOp, indices);

    // Slice off one corner to form a tri and a quad
    Pair<scalar> f;
    for (label i = 0; i < 2; ++ i)
    {
        f[i] = l[0]/(l[0] - l[i+1]);
    }
    if (l[0] > 0)
    {
        return triCutTri(a, p, f) + triCutQuad(b, p, f);
    }
    else
    {
        return triCutQuad(a, p, f) + triCutTri(b, p, f);
    }
}


template<class AboveOp, class BelowOp>
typename Foam::cut::opAddResult<AboveOp, BelowOp>::type Foam::triCut
(
    const FixedList<point, 3>& tri,
    const plane& p,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
)
{
    // Set the level set to the signed distance from the plane
    FixedList<scalar, 3> level;
    for (label i = 0; i < 3; ++ i)
    {
        level[i] = (tri[i] - p.refPoint()) & p.normal();
    }

    // Run the level set function
    return triCut(tri, level, aboveOp, belowOp);
}


template<class Point, class AboveOp, class BelowOp>
typename Foam::cut::opAddResult<AboveOp, BelowOp>::type Foam::tetCut
(
    const FixedList<Point, 4>& tet,
    const FixedList<scalar, 4>& level,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
)
{
    // Get the min and max over all four vertices and quick return if there is
    // no change of sign
    scalar levelMin = vGreat, levelMax = - vGreat;
    for (label i = 0; i < 4; ++ i)
    {
        levelMin = min(levelMin, level[i]);
        levelMax = max(levelMax, level[i]);
    }
    if (levelMin >= 0)
    {
        return aboveOp(tet) + belowOp();
    }
    if (levelMax <= 0)
    {
        return aboveOp() + belowOp(tet);
    }

    // Partition the level so that positive values are at the start. This is
    // like a single iteration of quick-sort, except that the pivot is a hard-
    // coded zero, rather than an element of the array. This can change the sign
    // of the tet.
    FixedList<label, 4> indices({0, 1, 2, 3});
    bool signChange = false;
    label i = 0, j = 3;
    while (true)
    {
        while (i < j && level[indices[i]] > 0)
        {
            i ++;
        }
        while (j > i && level[indices[j]] <= 0)
        {
            j --;
        }
        if (i == j)
        {
            break;
        }
        Swap(indices[i], indices[j]);
        signChange = !signChange;
    }

    // The number of vertices above the slice
    label n = i;

    // If there are more positives than negatives then reverse the order so that
    // the negatives are at the start
    if (n > 2)
    {
        n = 4 - n;
        for (label i = 0; i < 2; ++ i)
        {
            Swap(indices[i], indices[3-i]);
        }
    }

    // Correct the sign
    if (signChange)
    {
        Swap(indices[2], indices[3]);
    }

    // Permute the data
    const FixedList<Point, 4> p = tetReorder(tet, indices);
    const FixedList<scalar, 4> l = tetReorder(level, indices);
    AboveOp a = tetReorder(aboveOp, indices);
    BelowOp b = tetReorder(belowOp, indices);

    // Calculate the integrals above and below the level set
    if (n == 1)
    {
        // Slice off one corner to form a tet and a prism
        FixedList<scalar, 3> f;
        for (label i = 0; i < 3; ++ i)
        {
            f[i] = l[0]/(l[0] - l[i+1]);
        }
        if (l[0] > 0)
        {
            return tetCutTet(a, p, f) + tetCutPrism0(b, p, f);
        }
        else
        {
            return tetCutPrism0(a, p, f) + tetCutTet(b, p, f);
        }
    }
    else if (n == 2)
    {
        // Slice off two corners to form two prisms
        FixedList<scalar, 4> f;
        for (label i = 0; i < 2; ++ i)
        {
            for (label j = 0; j < 2; ++ j)
            {
                f[2*i+j] = l[i]/(l[i] - l[j+2]);
            }
        }
        if (l[0] > 0)
        {
            return tetCutPrism01(a, p, f) + tetCutPrism23(b, p, f);
        }
        else
        {
            return tetCutPrism23(a, p, f) + tetCutPrism01(b, p, f);
        }
    }

    FatalErrorInFunction
        << "The number of tet vertices above the level set should always be "
        << "either 1 or 2" << exit(FatalError);

    return aboveOp() + belowOp();
}


template<class AboveOp, class BelowOp>
typename Foam::cut::opAddResult<AboveOp, BelowOp>::type Foam::tetCut
(
    const FixedList<point, 4>& tet,
    const plane& p,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
)
{
    // Set the level set to the signed distance from the plane
    FixedList<scalar, 4> level;
    for (label i = 0; i < 4; ++ i)
    {
        level[i] = (tet[i] - p.refPoint()) & p.normal();
    }

    // Run the level set function
    return tetCut(tet, level, aboveOp, belowOp);
}

// ************************************************************************* //

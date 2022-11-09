/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "MPLICface.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MPLICface::MPLICface(const bool unweighted)
:
    unweighted_(unweighted)
{};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::MPLICface::cutFace
(
    const labelList& f,
    const labelList& faceEdges,
    const pointField& points,
    const boolList& isEdgeCutOld,
    boolList& isEdgeCut,
    label& faceEdgei,
    const UList<scalar>& pointsAlpha,
    const UList<vector>& pointsU,
    const label facei,
    const scalar target,
    const bool ow
)
{
    // Clear all the storage
    cutPoints_.clear();
    subPoints_.clear();
    cutEdges_.clear();
    subPointsU_.clear();

    // Direction of face circulators
    flipped_ = false;

    label fp;
    if (faceEdgei == -1)
    {
        const UIndirectList<scalar> fAlphas(pointsAlpha, f);

        // Starting from this point on the face
        fp = findMin(fAlphas);
        flipped_ = ow ? false : true;
    }
    else
    {
        // Finding first index of the edge in the face point list
        const label startPoint = faceEdgei;

        // Pick up starting point and decide direction of circulators
        const label fwd = f.fcIndex(startPoint);
        const label back = f.rcIndex(startPoint);
        const label lEdgeP = f[fwd] == f[f.fcIndex(faceEdgei)] ? fwd : back;

        // Starting from this point on the face
        fp = pointsAlpha[f[startPoint]] < pointsAlpha[f[lEdgeP]]
            ? startPoint : lEdgeP;

        if
        (
            !(
                (
                    lEdgeP == fwd
                 && pointsAlpha[f[startPoint]] < pointsAlpha[f[lEdgeP]]
                )
             || (
                    lEdgeP == back
                 && pointsAlpha[f[startPoint]] > pointsAlpha[f[lEdgeP]]
                )
             )
        )
        {
            flipped_ = true;
        }
    }

    label nextFp, edgei;
    if (flipped_)
    {
        nextFp = f.rcIndex(fp);
        edgei = nextFp;
    }
    else
    {
        nextFp = f.fcIndex(fp);
        edgei = fp;
    }

    forAll(f, i)
    {
        // Edge points field value
        const scalar& fl = pointsAlpha[f[fp]];
        const scalar& fk = pointsAlpha[f[nextFp]];

        const point& pL = points[f[fp]];
        const point& pK = points[f[nextFp]];

        const label edg = faceEdges[edgei];

        // Collect sub-point if the iso-value is bigger than target value
        if (fl >= target && cutPoints_.size() > 0)
        {
            subPoints_.append(pL);
            if (!unweighted_)
            {
                subPointsU_.append(pointsU[f[fp]]);
            }
        }

        // Cut the edge
        if
        (
            !isEdgeCutOld[edg]
         && ((fl < target && fk >= target) || (fl >= target && fk < target))
        )
        {
            const scalar coeff = (target - fk)/(fl - fk);
            const point cut = pK + coeff*(pL - pK);
            if (!unweighted_)
            {
                subPointsU_.append
                (
                    pointsU[f[nextFp]]
                  + coeff*(pointsU[f[fp]] - pointsU[f[nextFp]])
                );
            }
            cutPoints_.append(cut);
            subPoints_.append(cut);
            cutEdges_.append(edg);
            isEdgeCut[edg] =  true;
        }

        // Termination control
        if (cutPoints_.size() == 2)
        {
            faceEdgei = edgei;
            return 1;
        }

        if (flipped_)
        {
            fp = f.rcIndex(fp);
            nextFp = f.rcIndex(fp);
            edgei = nextFp;
        }
        else
        {
            fp = f.fcIndex(fp);
            nextFp = f.fcIndex(fp);
            edgei = fp;
        }
    }

    return 0;
}


Foam::label Foam::MPLICface::cutFace
(
    const UList<label>& f,
    const UList<point>& points,
    const UList<scalar>& pointsAlpha,
    const UList<vector>& pointsU,
    const scalar target,
    const bool ow
)
{
    // Clear all the storage
    cutPoints_.clear();
    subPoints_.clear();
    subPointsU_.clear();

    // Direction of face circulators
    flipped_ = ow ? false : true;

    const UIndirectList<scalar> fAlphas(pointsAlpha, f);

    // Starting from point with minimum alpha value
    label fp = findMin(fAlphas);

    // Second point index
    label nextFp = flipped_  ? f.rcIndex(fp) : f.fcIndex(fp);

    forAll(f, i)
    {
        // Edge points field value
        const scalar& fl = pointsAlpha[f[fp]];
        const scalar& fk = pointsAlpha[f[nextFp]];

        const point& pL = points[f[fp]];
        const point& pK = points[f[nextFp]];

        // Collect sub-point if the iso-value is bigger than target value
        if (fl >= target && cutPoints_.size() > 0)
        {
            subPoints_.append(pL);
            if (!unweighted_)
            {
                subPointsU_.append(pointsU[f[fp]]);
            }
        }

        // Cut the edge
        if ((fl < target && fk >= target) || (fl >= target && fk < target))
        {
            const scalar coeff = (target - fk)/(fl - fk);
            const point cut = pK + coeff*(pL - pK);
            if (!unweighted_)
            {
                subPointsU_.append
                (
                    pointsU[f[nextFp]]
                  + coeff*(pointsU[f[fp]] - pointsU[f[nextFp]])
                );
            }
            cutPoints_.append(cut);
            subPoints_.append(cut);
        }

        if (flipped_)
        {
            fp = f.rcIndex(fp);
            nextFp = f.rcIndex(fp);
        }
        else
        {
            fp = f.fcIndex(fp);
            nextFp = f.fcIndex(fp);
        }
    }

    // Simple cut
    if (cutPoints_.size() == 2)
    {
        return 1;
    }

    // Multiple cuts through the face
    else if (cutPoints_.size() > 2)
    {
        return -1;
    }

    // No cut
    else
    {
        return 0;
    }
}


// ************************************************************************* //

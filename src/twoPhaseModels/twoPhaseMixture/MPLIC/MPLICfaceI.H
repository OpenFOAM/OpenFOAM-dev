/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::MPLICface::alphaPhiU() const
{
    if (flipped_)
    {
        return -alphaPhiU(subPointsU_, subPoints_);
    }
    else
    {
        return alphaPhiU(subPointsU_, subPoints_);
    }
}


template<class VectorList, class PointList>
inline Foam::scalar Foam::MPLICface::alphaPhiU
(
    const VectorList& pointsU,
    const PointList& points
) const
{
    const point& baseP = points[0];
    const vector& baseU = pointsU[0];

    scalar alphaPhiU = 0;
    for(label i=1; i < points.size()-1; ++i)
    {
        // 1/2 Triangle area
        const vector area = (points[i] - baseP) ^ (points[i + 1] - baseP);

        // Calculate face area weighted velocity field
        // Average is missing 1/3 and area 1/2
        alphaPhiU += (baseU + pointsU[i] + pointsU[i + 1]) & area;
    }
    alphaPhiU /= 6.0;

    return alphaPhiU;
}


template<class VectorList, class PointList>
inline Foam::scalar Foam::MPLICface::alphaPhiU
(
    const VectorList& pointsU,
    const PointList& points,
    const labelList& f
) const
{
    const point& baseP = points[f[0]];
    const vector& baseU = pointsU[f[0]];

    scalar alphaPhiU = 0;
    for(label i=1; i < f.size()-1; ++i)
    {
        // 1/2 Triangle area
        const vector area = (points[f[i]] - baseP) ^ (points[f[i + 1]] - baseP);

        // Calculate face area weighted velocity field
        // Average is missing 1/3 and area 1/2
        alphaPhiU += (baseU + pointsU[f[i]] + pointsU[f[i + 1]]) & area;
    }
    alphaPhiU /= 6.0;

    return alphaPhiU;
}


inline const Foam::DynamicList<Foam::point>& Foam::MPLICface::cutPoints() const
{
    return cutPoints_;
}


inline const Foam::DynamicList<Foam::point>&
Foam::MPLICface::subPoints() const
{
    return subPoints_;
}


inline const Foam::DynamicList<Foam::label>& Foam::MPLICface::cutEdges() const
{
    return cutEdges_;
}


inline const Foam::vector Foam::MPLICface::Sf() const
{
    return face::area(subPoints_);
}


inline const Foam::vector Foam::MPLICface::Cf
(
    const vector& area
) const
{
    // If the face is a triangle, do a direct calculation
    if (subPoints_.size() == 3)
    {
        return (1.0/3.0)*(subPoints_[0] + subPoints_[1] + subPoints_[2]);
    }

    // For more complex faces, decompose into triangles ...

    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(subPoints_, pi)
    {
        pAvg += subPoints_[pi];
    }
    pAvg /= subPoints_.size();

    const vector sumAHat = normalised(area);

    // Compute the area-weighted sum of the triangle centres. Note use
    // the triangle area projected in the direction of the face normal
    // as the weight, *not* the triangle area magnitude. Only the
    // former makes the calculation independent of the initial estimate.
    scalar sumAn = 0;
    vector sumAnc = Zero;
    forAll(subPoints_, pi)
    {
        const point& p = subPoints_[pi];
        const point& pNext = subPoints_[subPoints_.fcIndex(pi)];

        const vector a = (pNext - p)^(pAvg - p);
        const vector c = p + pNext + pAvg;

        const scalar an = a & sumAHat;

        sumAn += an;
        sumAnc += an*c;
    }

    // Complete calculating centres and areas. If the face is too small
    // for the sums to be reliably divided then just set the centre to
    // the initial estimate.
    if (sumAn > vSmall)
    {
        return (1.0/3.0)*sumAnc/sumAn;
    }
    else
    {
        return pAvg;
    }
}


// ************************************************************************* //

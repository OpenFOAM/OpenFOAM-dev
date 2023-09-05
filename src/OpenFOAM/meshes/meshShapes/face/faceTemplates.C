/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "face.H"
#include "vectorAndError.H"

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

template<class PointField>
Foam::vector Foam::face::area(const PointField& ps)
{
    // If the face is a triangle, do a direct calculation
    if (ps.size() == 3)
    {
        return (1.0/2.0)*((ps[1] - ps[0])^(ps[2] - ps[0]));
    }

    // For more complex faces, decompose into triangles ...

    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(ps, pi)
    {
        pAvg += ps[pi];
    }
    pAvg /= ps.size();

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    vector sumA = Zero;
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vector a = (pNext - p)^(pAvg - p);

        sumA += a;
    }

    return (1.0/2.0)*sumA;
}


template<class PointField>
Foam::point Foam::face::centre(const PointField& ps)
{
    // The overhead associated with additionally calculating the area is small,
    // and the optimiser may even remove the additional expense all together,
    // so just use the face::areaAndCentre function and discard the area
    return areaAndCentre(ps).second();
}


template<class PointField>
Foam::Tuple2<Foam::vector, Foam::point> Foam::face::areaAndCentre
(
    const PointField& ps
)
{
    // If the face is a triangle, do a direct calculation
    if (ps.size() == 3)
    {
        return
            Tuple2<vector, point>
            (
                (1.0/2.0)*((ps[1] - ps[0])^(ps[2] - ps[0])),
                (1.0/3.0)*(ps[0] + ps[1] + ps[2])
            );
    }

    // For more complex faces, decompose into triangles ...

    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(ps, pi)
    {
        pAvg += ps[pi];
    }
    pAvg /= ps.size();

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    vector sumA = Zero;
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vector a = (pNext - p)^(pAvg - p);

        sumA += a;
    }
    const vector sumAHat = normalised(sumA);

    // Compute the area-weighted sum of the triangle centres. Note use
    // the triangle area projected in the direction of the face normal
    // as the weight, *not* the triangle area magnitude. Only the
    // former makes the calculation independent of the initial estimate.
    scalar sumAn = 0;
    vector sumAnc = Zero;
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vector a = (pNext - p)^(pAvg - p);
        const vector c = p + pNext + pAvg;

        const scalar an = a & sumAHat;

        sumAn += an;
        sumAnc += an*c;
    }

    // Complete calculating centres and areas. If the face is too small
    // for the sums to be reliably divided then just set the centre to
    // the point average.
    return
        Tuple2<vector, point>
        (
            (1.0/2.0)*sumA,
            sumAn > vSmall ? (1.0/3.0)*sumAnc/sumAn : pAvg
        );
}


template<class PointField>
Foam::Tuple2<Foam::vector, Foam::point> Foam::face::areaAndCentreStabilised
(
    const PointField& ps
)
{
    // If the face is a triangle then there are no stability concerns, so use
    // the un-stabilised algorithm
    if (ps.size() == 3)
    {
        return areaAndCentre(ps);
    }

    // As face::areaAndCentre
    point pAvg = Zero;
    forAll(ps, pi)
    {
        pAvg += ps[pi];
    }
    pAvg /= ps.size();

    // As face::areaAndCentre, but also track round-off error in the calculation
    vectorAndError sumA(Zero);
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vectorAndError a =
            (vectorAndError(pNext) - vectorAndError(p))
           ^(vectorAndError(pAvg) - vectorAndError(p));

        sumA += a;
    }
    const vectorAndError sumAHat = normalised(sumA);

    // As face::areaAndCentre, but also track round-off error in the calculation
    scalarAndError sumAn = 0;
    vector sumAnc(Zero);
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vectorAndError a =
            (vectorAndError(pNext) - vectorAndError(p))
           ^(vectorAndError(pAvg) - vectorAndError(p));
        const point c = p + pNext + pAvg;

        const scalarAndError an = a & sumAHat;

        sumAn += an;
        sumAnc += an.value*c;
    }

    // As face::areaAndCentre, but blend the calculated centroid and the point
    // average by the fraction of error in the calculated normal area
    Tuple2<vector, point> result((1.0/2.0)*sumA.value(), pAvg);
    if (sumAn.value >= vSmall)
    {
        const scalar f = min(max(sumAn.error/sumAn.value, 0), 1);
        result.second() = (1 - f)*(1.0/3.0)*sumAnc/sumAn.value + f*pAvg;
    }
    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::face::average
(
    const pointField& ps,
    const Field<Type>& fld
) const
{
    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return
            (1.0/3.0)
           *(
               fld[operator[](0)]
             + fld[operator[](1)]
             + fld[operator[](2)]
            );
    }

    // For more complex faces, decompose into triangles ...

    // Compute an estimate of the centre and the field average as the average
    // of the point values
    point pAvg = Zero;
    Type fldAvg = Zero;
    forAll(*this, pi)
    {
        pAvg += ps[operator[](pi)];
        fldAvg += fld[operator[](pi)];
    }
    pAvg /= size();
    fldAvg /= size();

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    vector sumA = Zero;
    forAll(*this, pi)
    {
        const point& p = ps[operator[](pi)];
        const point& pNext = ps[operator[](fcIndex(pi))];

        const vector a = (pNext - p)^(pAvg - p);

        sumA += a;
    }
    const vector sumAHat = normalised(sumA);

    // Compute the area-weighted sum of the average field values on each
    // triangle
    scalar sumAn = 0;
    Type sumAnf = Zero;
    forAll(*this, pi)
    {
        const point& p = ps[operator[](pi)];
        const point& pNext = ps[operator[](fcIndex(pi))];

        const vector a = (pNext - p)^(pAvg - p);
        const Type f =
            fld[operator[](pi)]
          + fld[operator[](fcIndex(pi))]
          + fldAvg;

        const scalar an = a & sumAHat;

        sumAn += an;
        sumAnf += an*f;
    }

    // Complete calculating the average. If the face is too small for the sums
    // to be reliably divided then just set the average to the initial
    // estimate.
    if (sumAn > vSmall)
    {
        return (1.0/3.0)*sumAnf/sumAn;
    }
    else
    {
        return fldAvg;
    }
}


// ************************************************************************* //

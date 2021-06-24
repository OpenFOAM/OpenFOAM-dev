/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "DynamicList.H"

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

template<class PointField>
Foam::point Foam::face::centre(const PointField& ps)
{
    // If the face is a triangle, do a direct calculation
    if (ps.size() == 3)
    {
        return (1.0/3.0)*(ps[0] + ps[1] + ps[2]);
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


template<class PointField>
Foam::vector Foam::face::area(const PointField& ps)
{
    // If the face is a triangle, do a direct calculation
    if (ps.size() == 3)
    {
        return 0.5*((ps[1] - ps[0])^(ps[2] - ps[0]));
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

    return 0.5*sumA;
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

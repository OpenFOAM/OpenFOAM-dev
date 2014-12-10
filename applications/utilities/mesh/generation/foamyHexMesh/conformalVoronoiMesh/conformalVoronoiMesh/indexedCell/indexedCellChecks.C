/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "plane.H"
#include "tetrahedron.H"
#include "pointConversion.H"
#include "CGALTriangulation3DKernel.H"


template<typename Cell>
Foam::scalar Foam::foamyHexMeshChecks::coplanarTet
(
    Cell& c,
    const scalar tol
)
{
    tetPointRef tet
    (
        topoint(c->vertex(0)->point()),
        topoint(c->vertex(1)->point()),
        topoint(c->vertex(2)->point()),
        topoint(c->vertex(3)->point())
    );

    const scalar quality = tet.quality();

    if (quality < tol)
    {
        return quality;
    }

    return 0;

//    plane triPlane
//    (
//        topoint(c->vertex(0)->point()),
//        topoint(c->vertex(1)->point()),
//        topoint(c->vertex(2)->point())
//    );
//
//    const scalar distance = triPlane.distance(topoint(c->vertex(3)->point()));
//
//    // Check if the four points are roughly coplanar. If they are then we
//    // cannot calculate the circumcentre. Better test might be the volume
//    // of the tet.
//    if (distance < tol)
//    {
//        return 0;
//    }
//
//    return distance;
}


template<typename Cell>
bool Foam::foamyHexMeshChecks::closePoints
(
    Cell& c,
    const scalar tol
)
{
    for (label v = 0; v < 4; ++v)
    {
        for (label vA = v + 1; vA < 4; ++vA)
        {
            if
            (
                mag
                (
                    topoint(c->vertex(v)->point())
                  - topoint(c->vertex(vA)->point())
                )
              < tol
            )
            {
                return true;
            }
        }
    }

    return false;
}


template<typename Cell>
bool Foam::foamyHexMeshChecks::smallVolume
(
    Cell& c,
    const scalar tol
)
{
    CGAL::Tetrahedron_3<baseK> tet
    (
        c->vertex(0)->point(),
        c->vertex(1)->point(),
        c->vertex(2)->point(),
        c->vertex(3)->point()
    );

    if (tet.volume() < tol)
    {
        return true;
    }

    return false;
}

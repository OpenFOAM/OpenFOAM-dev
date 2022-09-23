/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

Description
    Functions for constructing bounding spheres of lists of points

SourceFiles
    boundSphere.C
    boundSphereTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef boundSphere_H
#define boundSphere_H

#include "point.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Return whether or not the given sphere is valid
inline bool isValidBoundSphere(const Tuple2<point, scalar>& sphere);

//- Compute a sphere of four points or less where every point intersects the
//  sphere's surface
template<class PointField>
Tuple2<point, scalar> intersectBoundSphere
(
    const PointField& ps,
    const FixedList<label, 4>& pis,
    const label nPs
);

//- Compute a bounding sphere of four points or less
template<class PointField>
Tuple2<point, scalar> trivialBoundSphere
(
    const PointField& ps,
    const FixedList<label, 4>& pis,
    const label nPs
);

//- Compute a bounding sphere for an arbitrary number of points recursively
//  using Weizl's algorithm
template<class PointField>
Tuple2<point, scalar> weizlBoundSphere
(
    const PointField& ps,
    List<label>& pis,
    const label nPs,
    FixedList<label, 4>& boundaryPis,
    const label nBoundaryPs
);

//- Compute a bounding sphere for an arbitrary number of points, and given an
//  engine with which to randomise Weizl's algorithm
template<class PointField>
Tuple2<point, scalar> boundSphere(const PointField& ps, Random& rndGen);

//- Compute a bounding sphere for an arbitrary number of points
template<class PointField>
Tuple2<point, scalar> boundSphere(const PointField& ps);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "boundSphereTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

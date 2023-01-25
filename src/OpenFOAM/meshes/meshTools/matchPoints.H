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

Description
    Determine correspondence between points. See below.

SourceFiles
    matchPoints.C

\*---------------------------------------------------------------------------*/

#ifndef matchPoints_H
#define matchPoints_H

#include "scalarField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Function matchPoints Declaration
\*---------------------------------------------------------------------------*/

//- Determine correspondence between pointFields. Gets passed
//  local matching distance per point. Returns map from 0 to 1. Returns
//  true if found corresponding point in 1 for every point in 0; false
//  otherwise. Compares mag(pt - origin) and does proper comparison for points
//  with (almost) identical mag. pts0 can be subset of pts1.
bool matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<scalar>& matchDistance,
    const bool verbose,
    labelList& from0To1,
    const point& origin = point::zero
);


//- Supply pts0Dir and pts1Dir. They are directions associated with the points
//  e.g., a face normal associated with each face centre.
//  A match between similar points is only allowed if their directions are
//  equal and opposite
bool matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<point>& pts0Dir,
    const UList<point>& pts1Dir,
    const UList<scalar>& matchDistance,
    const bool verbose,
    labelList& from0To1,
    const point& origin = point::zero
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::pointToPointPlanarInterpolation

Description
    Interpolates between two sets of unstructured points using 2D Delaunay
    triangulation. Used in e.g. timeVaryingMapped bcs.

SourceFiles
    pointToPointPlanarInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef pointToPointPlanarInterpolation_H
#define pointToPointPlanarInterpolation_H

#include "FixedList.H"
#include "coordinateSystem.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pointToPointPlanarInterpolation Declaration
\*---------------------------------------------------------------------------*/

class pointToPointPlanarInterpolation
{
    // Private data

        //- Perturbation factor
        const scalar perturb_;

        //- Whether to use nearest point only (avoids triangulation, projection)
        const bool nearestOnly_;

        //- Coordinate system
        coordinateSystem referenceCS_;

        //- Number of source points (for checking)
        label nPoints_;

        //- Current interpolation addressing to face centres of underlying
        //  patch
        List<FixedList<label, 3>> nearestVertex_;

        //- Current interpolation factors to face centres of underlying
        //  patch
        List<FixedList<scalar, 3>> nearestVertexWeight_;

    // Private Member Functions

        //- Calculate a local coordinate system from set of points
        coordinateSystem calcCoordinateSystem(const pointField&) const;

        //- Calculate addressing and weights
        void calcWeights
        (
            const pointField& sourcePoints,
            const pointField& destPoints
        );

public:

    // Declare name of the class and its debug switch
    ClassName("pointToPointPlanarInterpolation");


    // Constructors

        //- Construct from 3D locations. Determines local coordinate system
        //  from sourcePoints and maps onto that. If nearestOnly skips any
        //  local coordinate system and triangulation and uses nearest vertex
        //  only
        pointToPointPlanarInterpolation
        (
            const pointField& sourcePoints,
            const pointField& destPoints,
            const scalar perturb,
            const bool nearestOnly = false
        );

        //- Construct from coordinate system and locations.
        pointToPointPlanarInterpolation
        (
            const coordinateSystem& referenceCS,
            const pointField& sourcePoints,
            const pointField& destPoints,
            const scalar perturb
        );


    // Member Functions

        //- Return the coordinateSystem
        const coordinateSystem& referenceCS() const
        {
            return referenceCS_;
        }

        //- Number of source points
        label sourceSize() const
        {
            return nPoints_;
        }

        //  patch
        const List<FixedList<label, 3>>& nearestVertex() const
        {
            return nearestVertex_;
        }

        //- Current interpolation factors to face centres of underlying
        //  patch
        const List<FixedList<scalar, 3>>& nearestVertexWeight() const
        {
            return nearestVertexWeight_;
        }

        //- Helper: extract words of times
        static wordList timeNames(const instantList&);

        //- Helper: find time. Return true if successful.
        static bool findTime
        (
            const instantList& times,
            const label startSampleTime,
            const scalar timeVal,
            label& lo,
            label& hi
        );

        //- Interpolate from field on source points to dest points
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& sourceFld) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointToPointPlanarInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

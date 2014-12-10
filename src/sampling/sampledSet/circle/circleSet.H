/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::circleSet

Description
    Samples along a circular path

SourceFiles
    circleSet.C

\*---------------------------------------------------------------------------*/

#ifndef circleSet_H
#define circleSet_H

#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class meshSearch;

/*---------------------------------------------------------------------------*\
                          Class circleSet Declaration
\*---------------------------------------------------------------------------*/

class circleSet
:
    public sampledSet
{
    // Private data

        // Circle definition

            //- Origin (x, y, z) in global cartesian co-ordinates
            point origin_;

            //- Axis of the circle
            vector circleAxis_;

            //- Point on circle (x, y, z) in global cartesian co-ordinates
            //  Defines start point
            point startPoint_;


        // Sampling definition

            //- Sampling interval in degrees about the origin
            scalar dTheta_;


    // Private Member Functions

        //- Samples all points in sampleCoords.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();

public:

    //- Runtime type information
    TypeName("circle");


    // Constructors

        //- Construct from components
        circleSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const point& origin,
            const vector& circleAxis,
            const point& startPoint,
            const scalar dTheta
        );

        //- Construct from dictionary
        circleSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    // Destructor

        virtual ~circleSet();


    // Member Functions

        //- Get reference point
        virtual point getRefPoint(const List<point>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

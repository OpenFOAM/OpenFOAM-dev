/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::sampledSets::lineCell

Description
    Cell-samples along a line at the mid-points in-between face-intersections

Usage
    \table
        Property    | Description                            | Req'd? | Default
        start       | The start point of the line            | yes    |
        end         | The end point of the line              | yes    |
        axis        | The coordinate axis that is written    | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type        lineCell;
        start       (0.5 0.6 0.5);
        end         (0.5 -0.3 -0.1);
        axis        x;
    }
    \endverbatim

SourceFiles
    lineCell.C

\*---------------------------------------------------------------------------*/

#ifndef lineCell_H
#define lineCell_H

#include "lineFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                         Class lineCell Declaration
\*---------------------------------------------------------------------------*/

class lineCell
:
    public sampledSet
{
    // Private data

        //- Start point
        const point start_;

        //- End point
        const point end_;


    // Private Member Functions

        //- Calculate all the sampling points
        virtual void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples and copies them into *this
        void genSamples();


public:

    //- Runtime type information
    TypeName("lineCell");


    // Static Member Functions

        //- Calculate the next mid point sample
        static void calcMidPointSample
        (
            const polyMesh& mesh,
            const point& prevPt,
            const label prevFace,
            const label prevSegment,
            const scalar prevCurveDist,
            const point& nextPt,
            const label nextFace,
            const label nextSegment,
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        );


    // Constructors

        //- Construct from dictionary
        lineCell
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~lineCell();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

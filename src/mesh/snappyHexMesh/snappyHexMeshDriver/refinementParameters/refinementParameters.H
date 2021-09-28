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

Class
    Foam::refinementParameters

Description
    Simple container to keep together refinement specific information.

SourceFiles
    refinementParameters.C

\*---------------------------------------------------------------------------*/

#ifndef refinementParameters_H
#define refinementParameters_H

#include "dictionary.H"
#include "labelList.H"
#include "point.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                    Class refinementParameters Declaration
\*---------------------------------------------------------------------------*/

class refinementParameters
{
public:

    // Public classes

        //- Class to hold the points to select cells inside and outside
        //  surface regions
        class cellSelectionPoints
        {
            //- Points inside the surface regions to selected cells
            const List<point> inside_;

            //- Points outside the surface region to deselect cells
            const List<point> outside_;


        public:

            //- Constructor
            cellSelectionPoints(const dictionary& dict);

            //- Return the points inside the surface regions to selected cells
            const List<point>& inside() const
            {
                return inside_;
            }

            //- Return the points outside the surface region to deselect cells
            const List<point>& outside() const
            {
                return outside_;
            }
        };


private:

    // Private Data

        //- Total number of cells
        const label maxGlobalCells_;

        //- Per processor max number of cells
        const label maxLocalCells_;

        //- When to stop refining
        const label minRefineCells_;

        //- Curvature
        scalar curvature_;

        //- Planarity criterion
        scalar planarAngle_;

        //- Number of layers between different refinement levels
        const label nBufferLayers_;

        //- Points to select cells inside and outside surface regions
        const cellSelectionPoints selectionPoints_;

        //- FaceZone faces allowed which have owner and neighbour in same
        //  cellZone?
        Switch allowFreeStandingZoneFaces_;

        //- Use old topology based problem-cell removal (cells with 8 points
        //  on surface)
        Switch useTopologicalSnapDetection_;

        //- Allowed load unbalance
        scalar maxLoadUnbalance_;

        Switch handleSnapProblems_;


public:

    // Constructors

        //- Construct from dictionary - new syntax
        refinementParameters(const dictionary& dict);

        //- Disallow default bitwise copy construction
        refinementParameters(const refinementParameters&) = delete;


    // Member Functions

        // Access

            //- Total number of cells
            label maxGlobalCells() const
            {
                return maxGlobalCells_;
            }

            //- Per processor max number of cells
            label maxLocalCells() const
            {
                return maxLocalCells_;
            }

            //- When to stop refining
            label minRefineCells() const
            {
                return minRefineCells_;
            }

            //- Curvature
            scalar curvature() const
            {
                return curvature_;
            }

            //- Angle when two intersections are considered to be planar
            scalar planarAngle() const
            {
                return planarAngle_;
            }

            //- Number of layers between different refinement levels
            label nBufferLayers() const
            {
                return nBufferLayers_;
            }

            //- Return the points to select cells inside and outside
            //  surface regions
            const cellSelectionPoints& selectionPoints() const
            {
                return selectionPoints_;
            }

            //- Are zone faces allowed only in between different cell zones
            //  or also just free standing?
            bool allowFreeStandingZoneFaces() const
            {
                return allowFreeStandingZoneFaces_;
            }

            //- Use old topology based problem-cell removal
            //  (cells with 8 points on surface)
            bool useTopologicalSnapDetection() const
            {
                return useTopologicalSnapDetection_;
            }

            //- Allowed load unbalance
            scalar maxLoadUnbalance() const
            {
                return maxLoadUnbalance_;
            }

            bool handleSnapProblems() const
            {
                return handleSnapProblems_;
            }


        // Other

            //- Checks that cells are in mesh. Returns cells they are in.
            labelList findCells(const polyMesh&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const refinementParameters&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

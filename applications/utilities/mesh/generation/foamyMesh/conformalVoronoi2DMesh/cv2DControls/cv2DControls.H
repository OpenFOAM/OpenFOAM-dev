/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::cv2DControls

Description
    Controls for the 2D CV mesh generator.

SourceFiles
    cv2DControls.C
    cv2DControlsI.H

\*---------------------------------------------------------------------------*/

#ifndef cv2DControls_H
#define cv2DControls_H

#include "Switch.H"
#include "dictionary.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class cv2DControls;

Ostream& operator<<(Ostream&, const cv2DControls&);


/*---------------------------------------------------------------------------*\
                        Class cv2DControls Declaration
\*---------------------------------------------------------------------------*/

class cv2DControls
{
    // Private Data

        //- Description of data_
        const dictionary& dict_;

        const dictionary& motionControl_;

        const dictionary& conformationControl_;


public:

        // Controls

            //- Minimum cell size below which protusions through the surface are
            //  not split
            scalar minCellSize_;

            //- Square of minCellSize
            scalar minCellSize2_;

            //- Maximum quadrant angle allowed at a concave corner before
            //  additional "mitering" lines are added
            scalar maxQuadAngle_;

            //- Near-wall region where cells are aligned with the wall
            scalar nearWallAlignedDist_;

            //- Square of nearWallAlignedDist
            scalar nearWallAlignedDist2_;

            //- Insert near-boundary point mirror or point-pairs
            Switch insertSurfaceNearestPointPairs_;

            //- Mirror near-boundary points rather than insert point-pairs
            Switch mirrorPoints_;

            //- Insert point-pairs vor dual-cell vertices very near the surface
            Switch insertSurfaceNearPointPairs_;

            Switch objOutput_;

            Switch meshedSurfaceOutput_;

            Switch randomiseInitialGrid_;

            scalar randomPerturbation_;

            label maxBoundaryConformingIter_;


        // Tolerances

            //- Maximum cartesian span of the geometry
            scalar span_;

            //- Square of span
            scalar span2_;

            //- Minimum edge-length of the cell size below which protusions
            //  through the surface are not split
            scalar minEdgeLen_;

            //- Square of minEdgeLen
            scalar minEdgeLen2_;

            //- Maximum notch size below which protusions into the surface are
            //  not filled
            scalar maxNotchLen_;

            //- Square of maxNotchLen
            scalar maxNotchLen2_;

            //- The minimum distance allowed between a dual-cell vertex
            // and the surface before a point-pair is introduced
            scalar minNearPointDist_;

            //- Square of minNearPoint
            scalar minNearPointDist2_;

            //- Distance between boundary conforming point-pairs
            scalar ppDist_;

            //- Square of ppDist
            scalar ppDist2_;


    // Constructors

        cv2DControls
        (
            const dictionary& controlDict,
            const boundBox& bb
        );

        //- Disallow default bitwise copy construction
        cv2DControls(const cv2DControls&) = delete;


    //- Destructor
    ~cv2DControls();


    // Member Functions

        // Access

            //- Return the minimum cell size
            inline scalar minCellSize() const;

            //- Return the square of the minimum cell size
            inline scalar minCellSize2() const;

            //- Return the maximum quadrant angle
            inline scalar maxQuadAngle() const;

            //- Return number of layers to align with the nearest wall
            inline scalar nearWallAlignedDist() const;

            //- Return square of nearWallAlignedDist
            inline scalar nearWallAlignedDist2() const;

            //- Return insertSurfaceNearestPointPairs Switch
            inline Switch insertSurfaceNearestPointPairs() const;

            //- Return mirrorPoints Switch
            inline Switch mirrorPoints() const;

            //- Return insertSurfaceNearPointPairs Switch
            inline Switch insertSurfaceNearPointPairs() const;

            //- Return the objOutput Switch
            inline Switch objOutput() const;

            //- Return the meshedSurfaceOutput Switch
            inline Switch meshedSurfaceOutput() const;

            //- Return the randomise initial point layout Switch
            inline Switch randomiseInitialGrid() const;

            //- Return the random perturbation factor
            inline scalar randomPerturbation() const;

            //- Return the maximum number of boundary conformation iterations
            inline label maxBoundaryConformingIter() const;

            //- Return the span
            inline scalar span() const;

            //- Return the span squared
            inline scalar span2() const;

            //- Return the minEdgeLen
            inline scalar minEdgeLen() const;

            //- Return the minEdgeLen squared
            inline scalar minEdgeLen2() const;

            //- Return the maxNotchLen
            inline scalar maxNotchLen() const;

            //- Return the maxNotchLen squared
            inline scalar maxNotchLen2() const;

            //- Return the minNearPointDist
            inline scalar minNearPointDist() const;

            //- Return the minNearPointDist squared
            inline scalar minNearPointDist2() const;

            //- Return the ppDist
            inline scalar ppDist() const;


        // Write

            //- Write controls to output stream.
            void write(Ostream& os) const;

            //- Ostream Operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const cv2DControls& s
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cv2DControls&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cv2DControlsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

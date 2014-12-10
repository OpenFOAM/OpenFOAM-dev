/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::wallPointYPlus

Description
    Holds information (coordinate and yStar) regarding nearest wall point.

    Used in VanDriest wall damping where the interest is in y+ but only
    needs to be calculated up to e.g. y+ < 200. In all other cells/faces
    the damping function becomes 1, since y gets initialized to GREAT and
    yStar to 1.

    Note: should feed the additional argument (yPlusCutoff) through as a
    template argument into FaceCellWave

SourceFiles
    wallPointYPlusI.H
    wallPointYPlus.C

\*---------------------------------------------------------------------------*/

#ifndef wallPointYPlus_H
#define wallPointYPlus_H

#include "wallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class wallPointYPlus Declaration
\*---------------------------------------------------------------------------*/

class wallPointYPlus
:
    public wallPointData<scalar>
{

private:

    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const wallPointYPlus& w2,
            const scalar tol,
            TrackingData& td
        );

public:

    // Static data members

        //- cut-off value for y+
        static scalar yPlusCutOff;


    // Constructors

        //- Construct null
        inline wallPointYPlus();

        //- Construct from origin, yStar, distance
        inline wallPointYPlus
        (
            const point& origin,
            const scalar yStar,
            const scalar distSqr
        );


    // Member Functions

        // Needed by FaceCellWave

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of cellI
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh& mesh,
                const label thisCellI,
                const label neighbourFaceI,
                const wallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of faceI
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const label neighbourCellI,
                const wallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of faceI
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const wallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );
};


//- Data associated with pointEdgePoint type as contiguous as underlying type
template<>
inline bool contiguous<wallPointYPlus>()
{
    return contiguous<wallPointData<scalar> >();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPointYPlusI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

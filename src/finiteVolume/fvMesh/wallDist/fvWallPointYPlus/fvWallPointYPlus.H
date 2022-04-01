/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::fvWallPointYPlus

Description
    Holds information (coordinate and yStar) regarding nearest wall point.

    Used in VanDriest wall damping where the interest is in y+ but only
    needs to be calculated up to e.g. y+ < 200. In all other cells/faces
    the damping function becomes 1, since y gets initialised to great and
    yStar to 1.

    Note: should feed the additional argument (yPlusCutoff) through as a
    template argument into FvFaceCellWave

SourceFiles
    fvWallPointYPlusI.H
    fvWallPointYPlus.C

\*---------------------------------------------------------------------------*/

#ifndef fvWallPointYPlus_H
#define fvWallPointYPlus_H

#include "fvWallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                       Class fvWallPointYPlus Declaration
\*---------------------------------------------------------------------------*/

class fvWallPointYPlus
:
    public fvWallPointData<scalar>
{
public:

    //- Class used to pass additional data in
    class trackData
    {
    public:

        //- Cut off distance
        scalar yPlusCutOff;
    };


private:

    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point& pt,
            const fvWallPointYPlus& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline fvWallPointYPlus();

        //- Construct from origin, yStar, distance
        inline fvWallPointYPlus
        (
            const point& origin,
            const scalar yStar,
            const scalar distSqr
        );


    // Member Functions

        // Needed by FvFaceCellWave

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of celli
            template<class TrackingData>
            inline bool updateCell
            (
                const fvMesh& mesh,
                const label thisCelli,
                const labelPair& neighbourPatchAndFacei,
                const fvWallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of facei
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const label neighbourCelli,
                const fvWallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of facei
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const fvWallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );
};


//- Data associated with pointEdgePoint type as contiguous as underlying type
template<>
inline bool contiguous<fvWallPointYPlus>()
{
    return contiguous<fvWallPointData<scalar>>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvWallPointYPlusI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

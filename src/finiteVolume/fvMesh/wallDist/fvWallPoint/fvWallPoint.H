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
    Foam::fvWallPoint

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    fvWallPointI.H
    fvWallPoint.C

\*---------------------------------------------------------------------------*/

#ifndef fvWallPoint_H
#define fvWallPoint_H

#include "wallPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvPatch;
class fvMesh;
class transformer;


/*---------------------------------------------------------------------------*\
                          Class fvWallPoint Declaration
\*---------------------------------------------------------------------------*/

class fvWallPoint
:
    public wallPoint
{
public:

    // Constructors

        //- Inherit constructors
        using wallPoint::wallPoint;


    // Member Functions

        // Needed by FvFaceCellWave

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const fvPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const fvMesh&,
                const label thisCelli,
                const labelPair& neighbourPatchAndFacei,
                const fvWallPoint& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh&,
                const labelPair& thisPatchAndFacei,
                const label neighbourCelli,
                const fvWallPoint& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh&,
                const labelPair& thisPatchAndFacei,
                const fvWallPoint& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );
};


//- Data associated with fvWallPoint type are contiguous
template<>
inline bool contiguous<fvWallPoint>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvWallPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

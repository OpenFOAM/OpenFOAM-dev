/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::WallInfo

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    WallInfoI.H

\*---------------------------------------------------------------------------*/

#ifndef WallInfo_H
#define WallInfo_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;

/*---------------------------------------------------------------------------*\
                          Class WallInfo Declaration
\*---------------------------------------------------------------------------*/

template<class WallLocation>
class WallInfo
:
    public WallLocation
{
public:

    // Constructors

        //- Inherit constructors
        using WallLocation::WallLocation;


    // Member Functions

        // Needed by FaceCellWave

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh& mesh,
                const WallInfo<WallLocation>& w2,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const WallInfo<WallLocation>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const WallInfo<WallLocation>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const WallInfo<WallLocation>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class wallPoint;

template<>
inline bool contiguous<WallInfo<wallPoint>>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class WallLocation, class Type>
class WallLocationData;

#define DefineContiguousWallLocationDataType(Type, nullArg)                   \
    template<>                                                                \
    inline bool contiguous<WallInfo<WallLocationData<wallPoint, Type>>>()     \
    {                                                                         \
        return true;                                                          \
    }

DefineContiguousWallLocationDataType(bool, );
DefineContiguousWallLocationDataType(label, );
FOR_ALL_FIELD_TYPES(DefineContiguousWallLocationDataType);

#undef DefineContiguousWallLocationDataType


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WallInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

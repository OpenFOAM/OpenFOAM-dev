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
    Foam::FvWallInfo

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    FvWallInfoI.H

\*---------------------------------------------------------------------------*/

#ifndef FvWallInfo_H
#define FvWallInfo_H

#include "wallPoint.H"
#include "wallFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvPatch;
class fvMesh;
class transformer;

/*---------------------------------------------------------------------------*\
                          Class FvWallInfoBase Declaration
\*---------------------------------------------------------------------------*/

template<class WallInfo, class Derived>
class FvWallInfoBase
:
    public WallInfo::template type<Derived>
{
public:

    // Constructors

        //- Inherit constructors
        using WallInfo::template type<Derived>::type;


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
                const FvWallInfoBase<WallInfo, Derived>& neighbourInfo,
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
                const FvWallInfoBase<WallInfo, Derived>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh&,
                const labelPair& thisPatchAndFacei,
                const FvWallInfoBase<WallInfo, Derived>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );
};

/*---------------------------------------------------------------------------*\
                          Class FvWallInfo Declaration
\*---------------------------------------------------------------------------*/

template<class WallInfo>
class FvWallInfo
:
    public FvWallInfoBase<WallInfo, FvWallInfo<WallInfo>>
{
    public:

        using
            FvWallInfoBase<WallInfo, FvWallInfo<WallInfo>>::
            FvWallInfoBase;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FvWallInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

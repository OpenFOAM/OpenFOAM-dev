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
    Foam::FvWallInfoYPlus

Description
    Holds information (coordinate and yStar) regarding nearest wall point.

    Used in VanDriest wall damping where the interest is in y+ but only
    needs to be calculated up to e.g. y+ < 200. In all other cells/faces
    the damping function becomes 1, since y gets initialised to great and
    yStar to 1.

SourceFiles
    FvWallInfoYPlusI.H

\*---------------------------------------------------------------------------*/

#ifndef FvWallInfoYPlus_H
#define FvWallInfoYPlus_H

#include "FvWallInfoData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class FvWallInfoYPlusBase Declaration
\*---------------------------------------------------------------------------*/

template<class WallInfo, class Derived>
class FvWallInfoYPlusBase
:
    public FvWallInfoDataBase<WallInfo, scalar, Derived>
{
public:

    //- Class used to pass additional data in
    class trackData
    {
    public:

        //- Cut off distance
        scalar yPlusCutOff;
    };


protected:

    friend class FvWallInfoBase<WallInfo, Derived>;


    // Protected Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point& pt,
            const FvWallInfoYPlusBase<WallInfo, Derived>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline FvWallInfoYPlusBase();

        //- Inherit constructors
        using FvWallInfoDataBase<WallInfo, scalar, Derived>::FvWallInfoDataBase;
};


/*---------------------------------------------------------------------------*\
                       Class FvWallInfoYPlus Declaration
\*---------------------------------------------------------------------------*/

template<class WallInfo>
class FvWallInfoYPlus
:
    public FvWallInfoYPlusBase<WallInfo, FvWallInfoYPlus<WallInfo>>
{
    public:

        using
            FvWallInfoYPlusBase<WallInfo, FvWallInfoYPlus<WallInfo>>::
            FvWallInfoYPlusBase;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FvWallInfoYPlusI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

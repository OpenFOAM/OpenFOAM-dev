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
    Foam::FvWallInfoData

Description
    Holds information (coordinate and normal) regarding nearest wall point.

    Is like FvWallInfo but transfer extra (passive) data.
    Used e.g. in wall distance calculation with wall reflection vectors.

SourceFiles
    FvWallInfoDataI.H

\*---------------------------------------------------------------------------*/

#ifndef FvWallInfoData_H
#define FvWallInfoData_H

#include "FvWallInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class WallInfo, class Type, class Derived>
class FvWallInfoDataBase;

template<class WallInfo, class Type, class Derived>
Istream& operator>>
(
    Istream&,
    FvWallInfoDataBase<WallInfo, Type, Derived>&
);

template<class WallInfo, class Type, class Derived>
Ostream& operator<<
(
    Ostream&,
    const FvWallInfoDataBase<WallInfo, Type, Derived>&
);


/*---------------------------------------------------------------------------*\
                        Class FvWallInfoDataBase Declaration
\*---------------------------------------------------------------------------*/

template<class WallInfo, class Type, class Derived>
class FvWallInfoDataBase
:
    public FvWallInfoBase<WallInfo, Derived>
{
    // Private Data

        //- Data at nearest wall center
        Type data_;


protected:

    // Protected Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const FvWallInfoDataBase<WallInfo, Type, Derived>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    friend class FvWallInfoBase<WallInfo, Derived>;


    typedef Type dataType;


    // Constructors

        //- Construct null
        inline FvWallInfoDataBase();

        //- Construct from data and other geometry
        template<class ... Geometry>
        inline FvWallInfoDataBase
        (
            const Type& data,
            const Geometry& ... geometry
        );


    // Member Functions

        // Access

            inline const Type& data() const;

            inline Type& data();

            template<class TrackingData>
            inline const Type& data(TrackingData& td) const;


        // Needed by meshWave

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const fvPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );


    // IOstream Operators

        friend Ostream& operator<< <WallInfo, Type, Derived>
        (
            Ostream&,
            const FvWallInfoDataBase<WallInfo, Type, Derived>&
        );

        friend Istream& operator>> <WallInfo, Type, Derived>
        (
            Istream&,
            FvWallInfoDataBase<WallInfo, Type, Derived>&
        );
};


/*---------------------------------------------------------------------------*\
                        Class FvWallInfoData Declaration
\*---------------------------------------------------------------------------*/

template<class WallInfo, class Type>
class FvWallInfoData
:
    public FvWallInfoDataBase<WallInfo, Type, FvWallInfoData<WallInfo, Type>>
{
    public:

        using
            FvWallInfoDataBase<WallInfo, Type, FvWallInfoData<WallInfo, Type>>::
            FvWallInfoDataBase;
};


#define DefineFvWallInfoType(Type, nullArg)                                    \
                                                                               \
    template<class WallInfo>                                                   \
    using CAT(FvWallInfo, CAPITALIZE(Type)) = FvWallInfoData<WallInfo, Type>;  \
                                                                               \
    template<>                                                                 \
    inline bool contiguous<CAT(FvWallInfo, CAPITALIZE(Type))<wallPoint>>()     \
    {                                                                          \
        return contiguous<FvWallInfo<wallPoint>>();                            \
    }

DefineFvWallInfoType(bool, nullArg)
DefineFvWallInfoType(label, nullArg)
FOR_ALL_FIELD_TYPES(DefineFvWallInfoType);

#undef DefineFvWallInfoType


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FvWallInfoDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

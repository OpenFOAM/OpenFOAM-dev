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
    Foam::PointEdgeLayerInfoData

Description
    Class to be used with PointEdgeWave which enumerates layers of points

SourceFiles
    PointEdgeLayerInfoDataI.H

\*---------------------------------------------------------------------------*/

#ifndef PointEdgeLayerInfoData_H
#define PointEdgeLayerInfoData_H

#include "pointEdgeLayerInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;

// Forward declaration of friend functions and operators
template<class Type>
class PointEdgeLayerInfoData;
template<class Type>
inline Ostream& operator<<(Ostream&, const PointEdgeLayerInfoData<Type>&);
template<class Type>
inline Istream& operator>>(Istream&, PointEdgeLayerInfoData<Type>&);

/*---------------------------------------------------------------------------*\
                          Class PointEdgeLayerInfoData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PointEdgeLayerInfoData
:
    public pointEdgeLayerInfo
{
    // Private Data

        //- Data
        Type data_;


public:

    // Constructors

        //- Construct null
        inline PointEdgeLayerInfoData();

        //- Construct given the point layer index
        inline PointEdgeLayerInfoData(const label pointLayer, const Type& data);


    // Member Functions

        // Access

            //- Return the data
            inline const Type& data() const;


        // Needed by PointEdgeWave

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchPointi,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of edge on point
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const label edgei,
                const PointEdgeLayerInfoData<Type>& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const PointEdgeLayerInfoData<Type>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No old information.
            template<class TrackingData>
            inline bool updatePoint
            (
                const PointEdgeLayerInfoData<Type>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgei,
                const label pointi,
                const PointEdgeLayerInfoData<Type>& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal
            (
                const PointEdgeLayerInfoData<Type>&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const PointEdgeLayerInfoData<Type>&) const;
        inline bool operator!=(const PointEdgeLayerInfoData<Type>&) const;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const PointEdgeLayerInfoData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            PointEdgeLayerInfoData<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DefineContiguousPointEdgeLayerInfoDataType(Type, nullArg)             \
    template<>                                                                \
    inline bool contiguous<PointEdgeLayerInfoData<Type>>()                    \
    {                                                                         \
        return true;                                                          \
    }

DefineContiguousPointEdgeLayerInfoDataType(bool, );
DefineContiguousPointEdgeLayerInfoDataType(label, );
FOR_ALL_FIELD_TYPES(DefineContiguousPointEdgeLayerInfoDataType);

#undef DefineContiguousPointEdgeLayerInfoDataType


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PointEdgeLayerInfoDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

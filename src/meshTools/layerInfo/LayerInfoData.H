/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::LayerInfoData

Description
    Class to be used with FaceCellWave which enumerates layers of cells and
    transports data through those layers

SourceFiles
    LayerInfoDataI.H

\*---------------------------------------------------------------------------*/

#ifndef LayerInfoData_H
#define LayerInfoData_H

#include "layerInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;

// Forward declaration of friend functions and operators
template<class Type>
class LayerInfoData;
template<class Type>
inline Ostream& operator<<(Ostream&, const LayerInfoData<Type>&);
template<class Type>
inline Istream& operator>>(Istream&, LayerInfoData<Type>&);

/*---------------------------------------------------------------------------*\
                        Class LayerInfoData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LayerInfoData
:
    public layerInfo
{
    // Private Data

        //- Data
        Type data_;


public:

    // Constructors

        //- Construct null
        inline LayerInfoData();

        //- Construct given the face layer index and direction
        inline LayerInfoData
        (
            const label faceLayer,
            const label direction,
            const Type& data
        );


    // Member Functions

        // Access

            //- Return the data
            inline const Type& data() const;


        // Needed by FaceCellWave

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label nbrFacei,
                const LayerInfoData<Type>& nbrInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label nbrCelli,
                const LayerInfoData<Type>& nbrInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const LayerInfoData<Type>& nbrInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test equality
            template<class TrackingData>
            inline bool equal
            (
                const LayerInfoData<Type>&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const LayerInfoData<Type>&) const;
        inline bool operator!=(const LayerInfoData<Type>&) const;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const LayerInfoData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            LayerInfoData<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DefineContiguousLayerInfoDataType(Type, nullArg)                      \
    template<>                                                                \
    inline bool contiguous<LayerInfoData<Type>>()                             \
    {                                                                         \
        return true;                                                          \
    }

DefineContiguousLayerInfoDataType(bool, );
DefineContiguousLayerInfoDataType(label, );
FOR_ALL_FIELD_TYPES(DefineContiguousLayerInfoDataType);

#undef DefineContiguousLayerInfoDataType


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LayerInfoDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

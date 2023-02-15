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
    Foam::PatchEdgeFacePointData

Description
    Transport of nearest point location, plus data, for use in
    PatchEdgeFaceWave

SourceFiles
    PatchEdgeFacePointDataI.H

\*---------------------------------------------------------------------------*/

#ifndef PatchEdgeFacePointData_H
#define PatchEdgeFacePointData_H

#include "patchEdgeFacePoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;

// Forward declaration of friend functions and operators
template<class Type>
class PatchEdgeFacePointData;
template<class Type>
Istream& operator>>(Istream&, PatchEdgeFacePointData<Type>&);
template<class Type>
Ostream& operator<<(Ostream&, const PatchEdgeFacePointData<Type>&);

/*---------------------------------------------------------------------------*\
                      Class PatchEdgeFacePointData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PatchEdgeFacePointData
:
    public patchEdgeFacePoint
{
    // Private Data

        //- Data at nearest wall centre
        Type data_;


public:

    // Constructors

        //- Construct null
        inline PatchEdgeFacePointData();

        //- Construct from data, origin, distance
        inline PatchEdgeFacePointData(const Type&, const point&, const scalar);


    // Member Functions

        // Access

            inline const Type& data() const;

            inline Type& data();

            template<class TrackingData>
            inline const Type& data(TrackingData& td) const;


        // Needed by meshWave

            //- Apply rotation matrix
            template<class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const label edgei,
                const label facei,
                const PatchEdgeFacePointData<Type>& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const PatchEdgeFacePointData<Type>& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const label facei,
                const label edgei,
                const PatchEdgeFacePointData<Type>& edgeInfo,
                const scalar tol,
                TrackingData& td
            );


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const PatchEdgeFacePointData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            PatchEdgeFacePointData<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DefineContiguousPatchEdgeFacePointDataType(Type, nullArg)             \
    template<>                                                                \
    inline bool contiguous<PatchEdgeFacePointData<Type>>()                    \
    {                                                                         \
        return true;                                                          \
    }

DefineContiguousPatchEdgeFacePointDataType(bool, );
DefineContiguousPatchEdgeFacePointDataType(label, );
FOR_ALL_FIELD_TYPES(DefineContiguousPatchEdgeFacePointDataType);

#undef DefineContiguousPatchEdgeFacePointDataType


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchEdgeFacePointDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

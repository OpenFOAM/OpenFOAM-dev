/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::patchEdgeFaceInfo

Description

SourceFiles
    patchEdgeFaceInfoI.H
    patchEdgeFaceInfo.C

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeFaceInfo_H
#define patchEdgeFaceInfo_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "tensor.H"
#include "pTraits.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class polyPatch;
class polyMesh;


// Forward declaration of friend functions and operators

class patchEdgeFaceInfo;

Istream& operator>>(Istream&, patchEdgeFaceInfo&);
Ostream& operator<<(Ostream&, const patchEdgeFaceInfo&);


/*---------------------------------------------------------------------------*\
                      Class patchEdgeFaceInfo Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeFaceInfo
{
    // Private Data

        //- Position of nearest wall center
        point origin_;

        //- Normal distance (squared) from point to origin
        scalar distSqr_;


    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const patchEdgeFaceInfo& w2,
            const scalar tol,
            TrackingData& td
        );

        //- Combine current with w2. Update distSqr, origin if w2 has smaller
        //  quantities and returns true.
        template<class TrackingData>
        inline bool update
        (
            const patchEdgeFaceInfo& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline patchEdgeFaceInfo();

        //- Construct from origin, distance
        inline patchEdgeFaceInfo(const point&, const scalar);


    // Member Functions

        // Access

            inline const point& origin() const;

            inline scalar distSqr() const;


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

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
                const label edgeI,
                const label facei,
                const patchEdgeFaceInfo& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const patchEdgeFaceInfo& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const label facei,
                const label edgeI,
                const patchEdgeFaceInfo& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const patchEdgeFaceInfo&, TrackingData& td) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const patchEdgeFaceInfo&) const;
        inline bool operator!=(const patchEdgeFaceInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchEdgeFaceInfo&);
        friend Istream& operator>>(Istream&, patchEdgeFaceInfo&);
};


//- Data associated with patchEdgeFaceInfo type are contiguous
template<>
inline bool contiguous<patchEdgeFaceInfo>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchEdgeFaceInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

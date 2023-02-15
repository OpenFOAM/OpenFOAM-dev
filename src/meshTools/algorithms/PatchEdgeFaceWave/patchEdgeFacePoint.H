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
    Foam::patchEdgeFacePoint

Description
    Transport of nearest point location for use in PatchEdgeFaceWave.

SourceFiles
    patchEdgeFacePointI.H

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeFacePoint_H
#define patchEdgeFacePoint_H

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

class patchEdgeFacePoint;

Istream& operator>>(Istream&, patchEdgeFacePoint&);
Ostream& operator<<(Ostream&, const patchEdgeFacePoint&);


/*---------------------------------------------------------------------------*\
                      Class patchEdgeFacePoint Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeFacePoint
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
            const patchEdgeFacePoint& w2,
            const scalar tol,
            TrackingData& td
        );

        //- Combine current with w2. Update distSqr, origin if w2 has smaller
        //  quantities and returns true.
        template<class TrackingData>
        inline bool update
        (
            const patchEdgeFacePoint& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline patchEdgeFacePoint();

        //- Construct from origin, distance
        inline patchEdgeFacePoint(const point&, const scalar);


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
                const label edgei,
                const label facei,
                const patchEdgeFacePoint& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const patchEdgeFacePoint& edgeInfo,
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
                const label edgei,
                const patchEdgeFacePoint& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal
            (
                const patchEdgeFacePoint&,
                TrackingData& td
            ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const patchEdgeFacePoint&) const;
        inline bool operator!=(const patchEdgeFacePoint&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const patchEdgeFacePoint&);
        inline friend Istream& operator>>(Istream&, patchEdgeFacePoint&);
};


//- Data associated with patchEdgeFacePoint type are contiguous
template<>
inline bool contiguous<patchEdgeFacePoint>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchEdgeFacePointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

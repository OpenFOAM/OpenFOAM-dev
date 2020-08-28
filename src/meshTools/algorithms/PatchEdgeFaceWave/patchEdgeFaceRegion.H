/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::patchEdgeFaceRegion

Description
    Transport of region for use in PatchEdgeFaceWave.

    Set element to -2 to denote blocked.

SourceFiles
    patchEdgeFaceRegionI.H
    patchEdgeFaceRegion.C

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeFaceRegion_H
#define patchEdgeFaceRegion_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "tensor.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;


// Forward declaration of friend functions and operators

class patchEdgeFaceRegion;

Istream& operator>>(Istream&, patchEdgeFaceRegion&);
Ostream& operator<<(Ostream&, const patchEdgeFaceRegion&);


/*---------------------------------------------------------------------------*\
                     Class patchEdgeFaceRegion Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeFaceRegion
{
    // Private Data

        //- region
        label region_;

    // Private Member Functions

        //- Combine current with w2. Update region_ if w2 has smaller
        //  quantities and returns true.
        template<class TrackingData>
        inline bool update
        (
            const patchEdgeFaceRegion& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline patchEdgeFaceRegion();

        //- Construct from origin, distance
        inline patchEdgeFaceRegion(const label);


    // Member Functions

        // Access

            inline label region() const;


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
                const indirectPrimitivePatch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const label edgeI,
                const label facei,
                const patchEdgeFaceRegion& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const patchEdgeFaceRegion& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const label facei,
                const label edgeI,
                const patchEdgeFaceRegion& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const patchEdgeFaceRegion&, TrackingData&) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const patchEdgeFaceRegion&) const;
        inline bool operator!=(const patchEdgeFaceRegion&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchEdgeFaceRegion&);
        friend Istream& operator>>(Istream&, patchEdgeFaceRegion&);
};


//- Data associated with patchEdgeFaceRegion type are contiguous
template<>
inline bool contiguous<patchEdgeFaceRegion>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchEdgeFaceRegionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

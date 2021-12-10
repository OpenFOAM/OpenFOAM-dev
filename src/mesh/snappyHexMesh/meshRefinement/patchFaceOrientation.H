/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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
    Foam::patchFaceOrientation

Description
    Transport of orientation for use in PatchEdgeFaceWave.

SourceFiles
    patchFaceOrientationI.H
    patchFaceOrientation.C

\*---------------------------------------------------------------------------*/

#ifndef patchFaceOrientation_H
#define patchFaceOrientation_H

#include "tensor.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;


// Forward declaration of friend functions and operators

class patchFaceOrientation;

Istream& operator>>(Istream&, patchFaceOrientation&);
Ostream& operator<<(Ostream&, const patchFaceOrientation&);


/*---------------------------------------------------------------------------*\
                    Class patchFaceOrientation Declaration
\*---------------------------------------------------------------------------*/

class patchFaceOrientation
{
    // Private Data

        //- Flip status
        label flipStatus_;

public:

    // Constructors

        //- Construct null
        inline patchFaceOrientation();

        //- Construct from components
        inline patchFaceOrientation(const label);


    // Member Functions

        //- Orientation
        inline label flipStatus() const;

        //- Reverse orientation
        inline void flip();


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
                const label edgei,
                const label facei,
                const patchFaceOrientation& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const patchFaceOrientation& edgeinfo,
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
                const label edgei,
                const patchFaceOrientation& edgeinfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const patchFaceOrientation&, TrackingData&) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const patchFaceOrientation&) const;
        inline bool operator!=(const patchFaceOrientation&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchFaceOrientation&);
        friend Istream& operator>>(Istream&, patchFaceOrientation&);
};


//- Data associated with patchFaceOrientation type are contiguous
template<>
inline bool contiguous<patchFaceOrientation>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchFaceOrientationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::patchEdgeFaceRegions

Description
    Transport of regions for use in PatchEdgeFaceWave.

    Set element to -1 to denote blocked.

SourceFiles
    patchEdgeFaceRegionsI.H
    patchEdgeFaceRegions.C

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeFaceRegions_H
#define patchEdgeFaceRegions_H

#include "labelList.H"
#include "scalar.H"
#include "tensor.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;


// Forward declaration of friend functions and operators

class patchEdgeFaceRegions;

Istream& operator>>(Istream&, patchEdgeFaceRegions&);
Ostream& operator<<(Ostream&, const patchEdgeFaceRegions&);


/*---------------------------------------------------------------------------*\
                           Class patchEdgeFaceRegions Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeFaceRegions
{
    // Private Data

        //- Region per point
        labelList regions_;

public:

    // Constructors

        //- Construct null
        inline patchEdgeFaceRegions();

        //- Construct from regions
        inline patchEdgeFaceRegions(const labelList&);

        //- Construct from regions (on edge)
        inline patchEdgeFaceRegions(const labelPair&);


    // Member Functions

        // Access

            inline const labelList& regions() const;


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Apply rotation matrix
            template<class Patch, class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const Patch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class Patch, class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const Patch& patch,
                const label edgeI,
                const label facei,
                const patchEdgeFaceRegions& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class Patch, class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const Patch& patch,
                const patchEdgeFaceRegions& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class Patch, class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const Patch& patch,
                const label facei,
                const label edgeI,
                const patchEdgeFaceRegions& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const patchEdgeFaceRegions&, TrackingData&) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const patchEdgeFaceRegions&) const;
        inline bool operator!=(const patchEdgeFaceRegions&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchEdgeFaceRegions&);
        friend Istream& operator>>(Istream&, patchEdgeFaceRegions&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchEdgeFaceRegionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

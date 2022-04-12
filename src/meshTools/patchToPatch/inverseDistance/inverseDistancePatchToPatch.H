/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::patchToPatches::inverseDistance

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    opposite face onto which its centre-normal ray projects, plus the immediate
    neighbours to that opposite face. The proportion of contribution from the
    different faces is calculated using inverse distance weighting.

SourceFiles
    inverseDistance.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDistancePatchToPatch_H
#define inverseDistancePatchToPatch_H

#include "patchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{

/*---------------------------------------------------------------------------*\
                          Class inverseDistance Declaration
\*---------------------------------------------------------------------------*/

class inverseDistance
:
    public patchToPatch
{
    // Private Member Data

        //- For each source face, the coupled target weights
        List<DynamicList<scalar>> srcWeights_;

        //- For each target face, the coupled source weights
        List<DynamicList<scalar>> tgtWeights_;


    // Private Member Functions

        //- Get the bound box for a source face
        virtual treeBoundBox srcBox
        (
            const face& srcFace,
            const pointField& srcPoints,
            const vectorField& srcPointNormals
        ) const;

        //- Return whether or not the face contains a point
        bool inside
        (
            const face& f,
            const pointField& ps,
            const point& p,
            const vector& r
        ) const;

        //- Intersect two faces
        virtual bool intersectFaces
        (
            const primitivePatch& patch,
            const primitivePatch& otherPatch,
            const label facei,
            const label otherFacei,
            DynamicList<label>& faceOtherFaces,
            DynamicList<scalar>& faceWeights
        ) const;

        //- Intersect two faces
        virtual bool intersectFaces
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch,
            const label srcFacei,
            const label tgtFacei
        );

        //- Initialisation
        virtual void initialise
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch
        );

        //- Send data that resulted from an intersection between the source
        //  patch and a distributed source-local-target patch back to the
        //  original target processes.
        virtual void rDistributeTgt
        (
            const primitiveOldTimePatch& tgtPatch,
            const distributionMap& tgtMap
        );

        //- Finalise the intersection
        virtual label finalise
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch,
            const transformer& tgtToSrc
        );


public:

    //- Runtime type information
    TypeName("inverseDistance");


    // Constructors

        //- Construct from components
        inverseDistance(const bool reverse);


    //- Destructor
    ~inverseDistance();


    // Member Functions

        //- For each source face, the coupled target weights
        virtual tmpNrc<List<DynamicList<scalar>>> srcWeights
        (
            const primitivePatch& srcPatch
        ) const;

        //- For each target face, the coupled source weights
        virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights
        (
            const primitivePatch& tgtPatch
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatches
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

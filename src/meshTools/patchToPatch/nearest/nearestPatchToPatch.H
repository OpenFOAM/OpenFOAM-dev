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
    Foam::patchToPatches::nearest

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    single nearest opposite face only.

SourceFiles
    nearest.C

\*---------------------------------------------------------------------------*/

#ifndef nearestPatchToPatch_H
#define nearestPatchToPatch_H

#include "patchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{

/*---------------------------------------------------------------------------*\
                          Class nearest Declaration
\*---------------------------------------------------------------------------*/

class nearest
:
    public patchToPatch
{
protected:

    // Private Member Data

        //- For each source face, the distance to its coupled target face
        List<scalar> srcDistances_;

        //- For each target face, the distance to its coupled source face
        List<scalar> tgtDistances_;


    // Private Member Functions

        //- Get the bound box for a source face
        virtual treeBoundBox srcBox
        (
            const face& srcFace,
            const pointField& srcPoints,
            const vectorField& srcPointNormals
        ) const;

        //- Intersect two faces
        bool intersectFaces
        (
            const primitivePatch& patch,
            const primitivePatch& otherPatch,
            const label facei,
            const label otherFacei,
            DynamicList<label>& faceOtherFaces,
            scalar& faceDistance
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

        //- Trim the local target patch so that only parts that actually
        //  intersect the source remain
        virtual labelList trimLocalTgt
        (
            const primitiveOldTimePatch& localTgtPatch
        );

        //- Send data that resulted from an intersection between the source
        //  patch and a distributed source-local-target patch back to the
        //  original target processes
        virtual void rDistributeTgt(const primitiveOldTimePatch& tgtPatch);

        //- For each source face, the coupled target weights
        virtual tmpNrc<List<DynamicList<scalar>>> srcWeights() const;

        //- For each target face, the coupled source weights
        virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights() const;


public:

    //- Runtime type information
    TypeName("nearest");


    // Constructors

        //- Construct from components
        nearest(const bool reverse);


    //- Destructor
    ~nearest();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatches
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

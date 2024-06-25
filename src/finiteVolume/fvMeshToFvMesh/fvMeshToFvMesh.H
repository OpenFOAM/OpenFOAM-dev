/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::fvMeshToFvMesh

Description

SourceFiles
    fvMeshToFvMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshToFvMesh_H
#define fvMeshToFvMesh_H

#include "meshToMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class fvMeshToFvMesh Declaration
\*---------------------------------------------------------------------------*/

class fvMeshToFvMesh
:
    public meshToMesh
{
private:

    // Private Typedefs

        //- Alias for surface boundary fields to reduce verbosity of method
        //  definitions below
        template<class Type>
        using SurfaceFieldBoundary =
            GeometricBoundaryField<Type, fvsPatchField, surfaceMesh>;


    // Private Data

        //- Reference to the source mesh
        const fvMesh& srcMesh_;

        //- Reference to the target mesh
        const fvMesh& tgtMesh_;


    // Private Member Functions

        //- Evaluate constraint types for the given vol field
        template<class Type>
        static void evaluateConstraintTypes(VolField<Type>& fld);


public:

    //- Run-time type information
    TypeName("fvMeshToFvMesh");


    // Constructors

        //- Construct from source and target meshes. If a patchMap is supplied,
        //  then interpolate between the specified patches. If not, then assume
        //  a consistent mesh with consistently named patches and interpolate
        //  1-to-1 between patches with the same name.
        fvMeshToFvMesh
        (
            const fvMesh& srcMesh,
            const fvMesh& tgtMesh,
            const word& engineType,
            const HashTable<word>& patchMap = NullObjectRef<HashTable<word>>()
        );


    //- Destructor
    virtual ~fvMeshToFvMesh();


    // Member Functions

        // Interpolation

            //- Interpolate a source vol field to the target with no left
            //  over values specified. If the interpolation weight sum is less
            //  than one for a face then they will be normalised. If the
            //  interpolation weight sum is zero for a face then that face's
            //  value will be NaN.
            template<class Type>
            tmp<VolField<Type>> srcToTgt
            (
                const VolField<Type>& srcFld
            ) const;

            //- Interpolate a source vol field to the target with left over
            //  values specified. If the interpolation weight sum is less than
            //  one for a face then the average will include the left over
            //  value multiplied by one minus the weight sum.
            template<class Type>
            tmp<VolField<Type>> srcToTgt
            (
                const VolField<Type>& srcFld,
                const VolField<Type>& leftOverTgtFld,
                const UList<wordRe>& tgtCuttingPatches
            ) const;

            //- Interpolate a source vol internal field to the target with no
            //  left over values specified. As the corresponding srcToTgt.
            template<class Type>
            tmp<VolInternalField<Type>> srcToTgt
            (
                const VolInternalField<Type>& srcFld
            ) const;

            //- Interpolate a source vol internal field to the target with left
            //  over values specified. As the corresponding srcToTgt.
            template<class Type>
            tmp<VolInternalField<Type>> srcToTgt
            (
                const VolInternalField<Type>& srcFld,
                const VolInternalField<Type>& leftOverTgtFld
            ) const;

            //- ...
            template<class Type>
            tmp<SurfaceFieldBoundary<Type>> srcToTgt
            (
                const SurfaceFieldBoundary<Type>& srcFld
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshToFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::fvMeshStitcher

Description
    Mesh manipulator that uses the intersection provided by the cyclic
    non-conformal poly patches to create non-conformal finite volume
    interfaces.

SourceFiles
    fvMeshStitcher.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshStitcher_H
#define fvMeshStitcher_H

#include "fvMesh.H"
#include "intersectionPatchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalCyclicFvPatch;

/*---------------------------------------------------------------------------*\
                        Class fvMeshStitcher Declaration
\*---------------------------------------------------------------------------*/

class fvMeshStitcher
{
protected:

    // Protected Typedefs

        //- Alias for surface boundary fields to reduce verbosity of method
        //  definitions below
        template<class Type>
        using SurfaceFieldBoundary =
            GeometricBoundaryField<Type, fvsPatchField, surfaceMesh>;

        //- Alias the patchToPatch intersection part struct
        using part = patchToPatches::intersection::part;

        //- Alias the patchToPatch intersection couple struct
        using couple = patchToPatches::intersection::couple;


private:

    // Private Static Data

        //- Prefix applied to the names of non-conformal fields which are
        //  stored for mapping purposes
        static const word nccFieldPrefix_;


    // Private Data

        //- Non-const fvMesh reference to allow update
        fvMesh& mesh_;


    // Private Member Functions

        // Intersection

            //- Perform intersections for the given non-conformal interface
            void intersectNonConformalCyclic
            (
                const nonConformalCyclicFvPatch& nccFvp,
                SurfaceFieldBoundary<label>& polyFacesBf,
                const tmp<SurfaceFieldBoundary<label>>& tOrigFacesNbrBf,
                const tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
                const tmp<SurfaceFieldBoundary<point>>& tOrigCfNbrBf,
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf,
                List<part>& origEdgeParts
            );

            //- Convert per-patch lists of edge parts into a single
            //  owner-orig-boundary list of edge parts that can be applied to
            //  the mesh
            List<part> calculateOwnerOrigBoundaryEdgeParts
            (
                const List<List<part>>& patchEdgeParts
            );

            //- Apply a list of owner-orig-boundary edge parts to
            //  non-orig-patch faces (both internal and boundary)
            void applyOwnerOrigBoundaryEdgeParts
            (
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf,
                const List<part>& ownerOrigBoundaryEdgeParts
            );

            //- Stabilise the faces that have had non-conformal coupled parts
            //  cut out from them
            void stabiliseOrigPatchFaces
            (
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf
            );

            //- Perform intersections for all non-conformal interfaces
            void intersectNonConformalCyclics
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                const tmp<SurfaceFieldBoundary<label>>& tOrigFacesNbrBf,
                const tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
                const tmp<SurfaceFieldBoundary<point>>& tOrigCfNbrBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            );

            //- Apply the given non-conformal topology (polyFacesBf) to the
            //  face geometry fields for the given patch type. Just set small
            //  stabilisation values. Don't actually do any cutting. Used for
            //  initialisation, or in situations when cutting cannot be done
            //  (e.g., when post processing a processor case).
            template<class NonConformalFvPatch>
            inline void createNonConformalStabilisationGeometry
            (
                const SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            );


        // Motion

            //- Return whether or not this stitcher supports a changing mesh
            virtual bool changing() const = 0;

            //- Initialise correction of the mesh fluxes
            virtual void initCorrectMeshPhi
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            ) = 0;

            //- Correct the mesh fluxes
            virtual void correctMeshPhi
            (
                const SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            ) = 0;


        // Field Mapping

            //- Store old times for all fields of a given type and class
            template<class Type, template<class> class GeoField>
            void storeOldTimeFields();

            //- Store old times for all fields of a given class
            template<template<class> class GeoField>
            void storeOldTimeFields();

            //- Resize the patch fields of a given type and class to match the
            //  mesh
            template<class Type, template<class> class GeoField>
            void resizePatchFields();

            //- Resize the patch fields of a given class to match the mesh
            template<template<class> class GeoField>
            void resizePatchFields();

            //- Pre-conform surface fields of a given type by separating NCC
            //  and original parts of non-conformal couplings and storing them
            //  in the database as separate fields
            template<class Type>
            void preConformSurfaceFields();

            //- Pre-conform surface fields by separating NCC and original
            //  parts of non-conformal couplings and storing them in the
            //  database as separate fields
            inline void preConformSurfaceFields();

            //- Post-non-conform surface fields of a given type by looking up
            //  NCC and original parts of non-conformal couplings and combining
            //  them into a single non-conformal boundary field
            template<class Type>
            void postNonConformSurfaceFields();

            //- Post-non-conform surface fields by looking up NCC and
            //  original parts of non-conformal couplings and combining them
            //  into a single non-conformal boundary field
            inline void postNonConformSurfaceFields();

            //- Evaluate all non-conformal vol patch fields of a given type
            template<class Type>
            void evaluateVolFields();

            //- Evaluate all non-conformal vol patch fields
            inline void evaluateVolFields();


protected:

    // Protected Member Functions

        // Field Mapping

            //- Reverse-map sum the values of a field
            template<class Type>
            static tmp<Field<Type>> fieldRMapSum
            (
                const Field<Type>& f,
                const label size,
                const labelUList& addr
            );

            //- Reverse-map sum the values of a field
            template<class Type>
            static tmp<Field<Type>> fieldRMapSum
            (
                const tmp<Field<Type>>& f,
                const label size,
                const labelUList& addr
            );

            //- Return the total non-conformal area associated with each
            //  original face
            inline tmp<SurfaceFieldBoundary<scalar>> getOrigNccMagSfb() const;

            //- Extract the non-conformal-coupled parts of the boundary field
            //  and store it on the conformal faces
            template<class Type>
            tmp<SurfaceFieldBoundary<Type>> conformalNccBoundaryField
            (
                const SurfaceFieldBoundary<Type>& fieldb
            ) const;

            //- Extract the original parts of the boundary field and store it
            //  on the conformal faces
            template<class Type>
            tmp<SurfaceFieldBoundary<Type>> conformalOrigBoundaryField
            (
                const SurfaceFieldBoundary<Type>& fieldb
            ) const;

            //- Combine non-conformal-coupled and original parts of the
            //  boundary field from the conformal faces to construct the
            //  complete non-conformal boundary field
            template<class Type>
            tmp<SurfaceFieldBoundary<Type>> nonConformalBoundaryField
            (
                const SurfaceFieldBoundary<Type>& nccFieldb,
                const SurfaceFieldBoundary<Type>& origFieldb
            ) const;

            //- Synchronise the boundary field by combining corresponding
            //  values across couplings with the given weightings
            template<class Type>
            tmp<SurfaceFieldBoundary<Type>> synchronisedBoundaryField
            (
                const SurfaceFieldBoundary<Type>& fieldb,
                const bool flip,
                const scalar ownerWeight,
                const scalar neighbourWeight
            ) const;

            //- Synchronise the boundary field by combining corresponding
            //  values across couplings with equal weightings
            template<class Type>
            tmp<SurfaceFieldBoundary<Type>> synchronisedBoundaryField
            (
                const SurfaceFieldBoundary<Type>& fieldb
            ) const;


        // Checking

            //- Return the non-dimensional cell openness for debugging/checking
            tmp<DimensionedField<scalar, volMesh>> openness() const;

            //- Return the non-dimensional cell openness for debugging/checking
            tmp<DimensionedField<scalar, volMesh>> volumeConservationError
            (
                const surfaceScalarField& meshPhi
            ) const;

            //- Return the non-dimensional cell openness for debugging/checking
            tmp<DimensionedField<scalar, volMesh>> volumeConservationError
            (
            ) const;


public:

    //- Runtime type information
    TypeName("fvMeshStitcher");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvMeshStitcher,
            fvMesh,
            (fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        explicit fvMeshStitcher(fvMesh&);

        //- Disallow default bitwise copy construction
        fvMeshStitcher(const fvMeshStitcher&) = delete;


    // Selectors

        //- Select, construct and return the fvMeshStitcher
        static autoPtr<fvMeshStitcher> New(fvMesh&, const bool changing);


    //- Destructor
    virtual ~fvMeshStitcher();


    // Member Functions

        //- Return the fvMesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Does this stitcher do anything?
        bool stitches() const;

        //- Update local data for topology changes
        virtual void updateMesh(const polyTopoChangeMap&);

        //- Update local data for mesh motion
        virtual void movePoints();

        //- Disconnect the mesh by removing faces from the nonConformalCyclics
        virtual bool disconnect
        (
            const bool changing,
            const bool geometric
        );

        //- Connect the mesh by adding faces into the nonConformalCyclics
        virtual bool connect
        (
            const bool changing,
            const bool geometric,
            const bool load
        );

        //- Write the mover state
        virtual bool write(const bool write = true) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshStitcher&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshStitcherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

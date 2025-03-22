/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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
#include "HashPtrTable.H"
#include "intersectionPatchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalCyclicFvPatch;
class nonConformalMappedWallFvPatch;

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
        using SurfaceFieldBoundary = GeometricBoundaryField<Type, surfaceMesh>;

        //- Alias the patchToPatch intersection part struct
        using part = patchToPatches::intersection::part;

        //- Alias the patchToPatch intersection couple struct
        using couple = patchToPatches::intersection::couple;


private:

    // Private Static Data

        //- Minimum projected volume fraction that results in a warning
        static const scalar minWarnProjectedVolumeFraction_;


    // Private Data

        //- Non-const fvMesh reference to allow update
        fvMesh& mesh_;

        //- Flag to indicate whether this region is ready to be stitched
        bool ready_;

        //- Cache of loaded region poly face IO objects
        HashPtrTable<IOobject> regionPolyFacesBfIOs_;

        //- Cache of loaded region poly face boundary fields
        HashPtrTable<SurfaceFieldBoundary<label>> regionPolyFacesBfs_;


    // Private Member Functions

        // Regions

            //- Construct a table of connected region names. Recurses to
            //  find neighbours-of-neighbours and so on. Includes this region.
            static void regionNames(const fvMesh& mesh, wordHashSet& names);

            //- Get a list of the connected region names
            wordList regionNames() const;

            //- Get a list of the connected region meshes
            UPtrList<const fvMesh> regionMeshes() const;

            //- Get a list of the connected region meshes
            UPtrList<fvMesh> regionMeshes();


        // Topology access

            //- Get the IO for the poly faces boundary field for a given mesh
            static IOobject polyFacesBfIO(const fvMesh& mesh);

            //- Load the poly faces boundary field for this mesh, and also make
            //  sure that all neighbouring regions' poly faces boundary fields
            //  are available either by lookup or in the region cache
            bool loadPolyFacesBf
            (
                IOobject& polyFacesBfIO,
                SurfaceFieldBoundary<label>& polyFacesBf
            );

            //- Access the poly faces boundary field for a region
            const SurfaceFieldBoundary<label>& getPolyFacesBf
            (
                const word& regionName
            ) const;

            //- Get fields of neighbouring original face properties
            void getOrigNbrBfs
            (
                const SurfaceFieldBoundary<label>& polyFacesBf,
                const SurfaceFieldBoundary<vector>& SfBf,
                const SurfaceFieldBoundary<vector>& CfBf,
                tmp<SurfaceFieldBoundary<label>>& tOrigFacesNbrBf,
                tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
                tmp<SurfaceFieldBoundary<point>>& tOrigCfNbrBf
            ) const;


        // Intersection

            //- Unpack the patchToPatch addressing into lists of indices (fixed
            //  lists of 3 labels; owner face, neighbour face, couple index).
            //  These will be used to create the non-conformal faces, so sort
            //  them to make sure the non-conformal interfaces are ordered.
            static List<List<FixedList<label, 3>>> procFacesToIndices
            (
                const List<List<remote>>& faceOtherProcFaces,
                const bool owner
            );

            //- If addressing has been provided, then modify the indices to
            //  match. When a coupling has to be added, the couple index is set
            //  to -1. This indicates that there is no geometry in the
            //  patchToPatch engine for this coupling, and for a small
            //  stabilisation value to be used instead.
            static void matchIndices
            (
                const SurfaceFieldBoundary<label>& polyFacesBf,
                const SurfaceFieldBoundary<label>& origFacesNbrBf,
                List<List<FixedList<label, 3>>>& indices,
                const fvPatch& ncFvp,
                const polyPatch& origPp,
                const labelList& patchis,
                const labelList& patchOffsets,
                const labelList& patchSizes,
                const bool owner
            );

            //- Return the number of valid indices in the given list
            static label nValidIndices
            (
                const List<FixedList<label, 3>>& indices
            );

            //- Create couplings by transferring geometry from the original to
            //  the non-conformal patches
            static void createCouplings
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf,
                const tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
                const tmp<SurfaceFieldBoundary<vector>>& tOrigCfNbrBf,
                const List<List<FixedList<label, 3>>>& indices,
                const List<DynamicList<couple>>& couples,
                const polyPatch& origPp,
                const labelList& patchis,
                const labelList& patchOffsets,
                const bool owner
            );

            //- Add error geometry to the original patches and store the edge
            //  parts
            static void createErrorAndEdgeParts
            (
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf,
                List<part>& origEdgeParts,
                const patchToPatches::intersection& intersection,
                const polyPatch& origPp
            );

            //- Perform intersections for the given non-conformal cyclic
            void intersectNonConformalCyclic
            (
                const nonConformalCyclicFvPatch& nccFvp,
                SurfaceFieldBoundary<label>& polyFacesBf,
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf,
                const tmp<SurfaceFieldBoundary<label>>& tOrigFacesNbrBf,
                const tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
                const tmp<SurfaceFieldBoundary<point>>& tOrigCfNbrBf,
                List<part>& origEdgeParts
            ) const;

            //- Perform intersections for the given non-conformal mapped
            void intersectNonConformalMappedWall
            (
                const nonConformalMappedWallFvPatch& ncmwFvp,
                SurfaceFieldBoundary<label>& polyFacesBf,
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf,
                const tmp<SurfaceFieldBoundary<label>>& tOrigFacesNbrBf,
                const tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
                const tmp<SurfaceFieldBoundary<point>>& tOrigCfNbrBf,
                List<part>& origEdgeParts
            ) const;

            //- Convert per-patch lists of edge parts into a single
            //  owner-orig-boundary list of edge parts that can be applied to
            //  the mesh
            List<part> calculateOwnerOrigBoundaryEdgeParts
            (
                const List<List<part>>& patchEdgeParts
            ) const;

            //- Apply a list of owner-orig-boundary edge parts to
            //  non-orig-patch faces (both internal and boundary)
            void applyOwnerOrigBoundaryEdgeParts
            (
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf,
                const List<part>& ownerOrigBoundaryEdgeParts
            ) const;

            //- Stabilise the faces that have had non-conformal coupled parts
            //  cut out from them
            void stabiliseOrigPatchFaces
            (
                SurfaceFieldBoundary<vector>& SfBf,
                SurfaceFieldBoundary<vector>& CfBf
            ) const;

            //- Perform intersections for all non-conformal interfaces. For
            //  patches that can be coupled, this creates the intersection
            //  geometry. For patches that cannot, this just sets small
            //  stabilisation values for the given topology (if any).
            void intersect
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf,
                const boolList& patchCoupleds,
                const bool matchTopology
            ) const;


        // Motion

            //- Return whether or not this stitcher supports a changing mesh
            virtual bool changing() const = 0;

            //- Correct the mesh fluxes following a conform operation
            virtual void conformCorrectMeshPhi(surfaceScalarField& phi) = 0;

            //- Initialise correction of the mesh fluxes
            virtual void createNonConformalCorrectMeshPhiGeometry
            (
                SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf
            ) = 0;

            //- Correct the mesh fluxes following an unconform operation
            virtual void unconformCorrectMeshPhi
            (
                const SurfaceFieldBoundary<label>& polyFacesBf,
                surfaceVectorField& SfSf,
                surfaceVectorField& CfSf,
                surfaceScalarField& phi
            ) = 0;


        // Connection

            //- Disconnect the mesh by removing faces from the
            //  nonConformalCyclics
            bool disconnectThis
            (
                const bool changing,
                const bool geometric
            );

            //- Connect the mesh by adding faces into the nonConformalCyclics
            bool connectThis
            (
                const bool changing,
                const bool geometric,
                const bool load
            );


        // Field Mapping

            //- Resize the patch fields of a given type and class to match the
            //  sizes of the patches in the mesh
            template<class Type, template<class> class GeoField>
            void resizePatchFields();

            //- As above, for all types
            template<template<class> class GeoField>
            void resizePatchFields();

            //- Pre-conform surface fields of a given type. This maps all the
            //  non-conformal values to the original patch faces so that they
            //  can be mapped across subsequent mesh changes. See
            //  conformedFvsPatchField::conform for details.
            template<class Type>
            void preConformSurfaceFields();

            //- Calls pre-conform surface field functions for all types
            void preConformSurfaceFields();

            //- Pre-conform the volume fields of a given type. This maps all the
            //  non-conformal values to the original patch faces so that they
            //  can be mapped across subsequent mesh changes. See
            //  conformedFvPatchField::conform for details.
            template<class Type>
            void preConformVolFields();

            //- Calls pre-conform volume field functions for all types
            void preConformVolFields();

            //- Post-unconform surface fields of a given type. This takes all
            //  the non-conformal values which were mapped onto the conformal
            //  faces by the pre-conform call and maps them back into the new
            //  non-conformal faces. See conformedFvsPatchField::unconform for
            //  details.
            template<class Type>
            void postUnconformSurfaceFields();

            //- Calls post-unconform surface field functions for all types
            void postUnconformSurfaceFields();

            //- Post-unconform volume fields of a given type. This takes all
            //  the non-conformal values which were mapped onto the conformal
            //  faces by the pre-conform call and maps them back into the new
            //  non-conformal faces. See conformedFvPatchField::unconform for
            //  details.
            template<class Type>
            void postUnconformVolFields();

            //- Evaluate all non-conformal vol patch fields of a given type
            template<class Type>
            void postUnconformEvaluateVolFields();

            //- Calls post-unconform volume field functions for all types
            void postUnconformVolFields();


protected:

    // Protected Member Functions

        // Checking

            //- Determine which patches are coupled; i.e., for which
            //  non-conformal patches are the remote or neighbouring meshes
            //  available
            boolList patchCoupleds() const;

            //- Is the connection "geometric", or has the topology just been
            //  loaded and stabilisation geometry applied?
            bool geometric() const;

            //- Return the non-dimensional cell openness for debugging/checking
            tmp<DimensionedField<scalar, volMesh>> openness() const;

            //- Return the non-dimensional old-time volume conservation error
            //  for a specified old-time index for debugging/checking
            tmp<DimensionedField<scalar, volMesh>>
                volumeConservationError(const label n) const;

            //- Return the projected volume fraction for debugging/checking
            tmp<DimensionedField<scalar, volMesh>>
                projectedVolumeFraction() const;


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
        static autoPtr<fvMeshStitcher> New(fvMesh&, bool changing);


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

        //- Disconnect the mesh by removing faces from the nonConformalCyclics
        bool disconnect
        (
            const bool changing,
            const bool geometric
        );

        //- Connect the mesh by adding faces into the nonConformalCyclics
        bool connect
        (
            const bool changing,
            const bool geometric,
            const bool load
        );

        //- Re-compute the connection. Topology is preserved. Permits a change
        //  in whether or not the connection is "geometric".
        void reconnect(const bool geometric) const;

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Write the stitcher state
        virtual bool write(const bool write = true) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshStitcher&) = delete;
};


// Template Specialisations

template<>
void fvMeshStitcher::postUnconformSurfaceFields<Foam::vector>();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshStitcherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

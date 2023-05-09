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
    Foam::mappedPatchBase

Description
    Engine which provides mapping between two patches

    Example:
    \verbatim
        // The name of the region to map from. Optional. Defaults to the same
        // region as the patch.
        neighbourRegion region0;

        // The name of the patch to map from
        neighbourPatch movingWall;

        // Couple group to specify the region and patch to map from. This is an
        // alternative to specifying neighbourRegion and neighbourPatch
        // directly, as shown above.
        coupleGroup baffleGroup;
    \endverbatim

SourceFiles
    mappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchBase_H
#define mappedPatchBase_H

#include "pointField.H"
#include "patchToPatch.H"
#include "coupleGroupIdentifier.H"
#include "cyclicTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;
class distributionMap;

/*---------------------------------------------------------------------------*\
                       Class mappedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchBase
{
protected:

    // Protected data

        //- Patch to map to
        const polyPatch& patch_;

        //- Couple group for the region/patch to map from
        const coupleGroupIdentifier coupleGroup_;

        //- Name of the region to map from
        mutable word nbrRegionName_;

        //- Name of the patch to map from
        mutable word nbrPatchName_;

        //- The transformation between the patches
        mutable cyclicTransform transform_;

        //- Are we using the tree mapping method, or a patch-to-patch
        //  intersection engine? The former is typically faster on small,
        //  pairs of patches with identical meshes. The latter is better
        //  parallelised and provides options for connecting patches with
        //  differing surface meshes.
        const bool usingTree_;

        //- Distributor (if using tree)
        mutable autoPtr<distributionMap> treeMapPtr_;

        //- Indices of the neighbouring patch faces who's values
        //  have to be supplied to the distribution map (if using tree)
        mutable labelList treeNbrPatchFaceIndices_;

        //- Is the patch-to-patch intersection engine up to date? (if using
        //  patch-to-patch)
        mutable bool patchToPatchIsValid_;

        //- Patch-to-patch intersection engine (if using patch-to-patch)
        mutable autoPtr<patchToPatch> patchToPatchPtr_;

        //- Default Matching tolerance
        static const scalar defaultMatchTol_;

        //- Matching tolerance
        const scalar matchTol_;

        //- Do we need to re-calculate the mapping if mesh motion takes place?
        //  Defaults to true.
        const bool reMapAfterMove_;

        //- Do we need to re-map the neighbour because this patch moved?
        mutable bool reMapNbr_;


    // Protected Member Functions

        //- Get the face-areas for this patch
        virtual tmp<vectorField> patchFaceAreas() const;

        //- Get the face-centres for this patch
        virtual tmp<pointField> patchFaceCentres() const;

        //- Get the local points for this patch
        virtual tmp<pointField> patchLocalPoints() const;

        //- Get the face-areas for the neighbour patch
        tmp<vectorField> nbrPatchFaceAreas() const;

        //- Get the face-centres for the neighbour patch
        tmp<pointField> nbrPatchFaceCentres() const;

        //- Get the local points for the neighbour patch
        tmp<pointField> nbrPatchLocalPoints() const;

        //- Calculate mapping
        void calcMapping() const;

        //- Return whether or not the neighbour patch is of mapped type
        inline bool nbrPatchIsMapped() const;

        //- Get the mapped neighbouring patch
        inline const mappedPatchBase& nbrMappedPatch() const;

        //- Determine whether the neighbouring patch maps from this patch in
        //  the same way as this patch maps from its neighbour
        inline bool symmetric() const;


public:

    //- Runtime type information
    TypeName("mappedPatchBase");


    // Constructors

        //- Construct from patch
        mappedPatchBase(const polyPatch&);

        //- Construct from components
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& nbrRegionName,
            const word& nbrPatchName,
            const cyclicTransform& transform
        );

        //- Construct from dictionary
        mappedPatchBase
        (
            const polyPatch& pp,
            const dictionary& dict,
            const bool defaultTransformIsNone
        );

        //- Construct as copy, resetting patch
        mappedPatchBase(const polyPatch&, const mappedPatchBase&);


    //- Destructor
    virtual ~mappedPatchBase();


    // Member Functions

        // Access

            //- Name of the region to map from
            inline const word& nbrRegionName() const;

            //- Is the neighbour region the same as for this patch?
            inline bool sameRegion() const;

            //- Is the neighbour patch the same as this patch?
            inline bool samePatch() const;

            //- Is the neighbour patch the same as this patch with no transform?
            inline bool sameUntransformedPatch() const;

            //- Name of the patch to map from
            inline const word& nbrPatchName() const;

            //- Get the mesh for the region to map from
            const polyMesh& nbrMesh() const;

            //- Get the patch to map from
            const polyPatch& nbrPolyPatch() const;


        // Casting

            //- Cast the given polyPatch to a mappedPatchBase. Handle errors.
            static const mappedPatchBase& getMap(const polyPatch& patch);

            //- Restrict use of the mapper to certain configurations
            struct from
            {
                static const label any = 0;
                static const label sameRegion = 1;
                static const label differentRegion = 2;
                static const label differentPatch = 4;
            };

            //- Validate that the map exists and is appropriate for the given
            //  set of permitted configurations
            template<class PatchFieldType, class FieldType>
            static void validateMapForField
            (
                const PatchFieldType& field,
                const FieldType& iF,
                const dictionary& context,
                const label froms = from::any
            );

            //- Validate that the map is appropriate for the given
            //  set of permitted configurations
            template<class PatchFieldType, class FieldType>
            void validateForField
            (
                const PatchFieldType& field,
                const FieldType& iF,
                const dictionary& context,
                const label froms = from::any
            ) const;


        // Edit

            //- Clear out data on mesh change
            virtual void clearOut();


        // Distribute

            //- Map/interpolate the neighbour patch field to this patch
            template<class Type>
            tmp<Field<Type>> fromNeighbour(const Field<Type>& nbrFld) const;

            //- Map/interpolate the neighbour patch field to this patch
            template<class Type>
            tmp<Field<Type>> fromNeighbour
            (
                const tmp<Field<Type>>& nbrFld
            ) const;

            //- Map/interpolate patch field to the neighbour patch
            template<class Type>
            tmp<Field<Type>> toNeighbour(const Field<Type>& fld) const;

            //- Map/interpolate patch field to the neighbour patch
            template<class Type>
            tmp<Field<Type>> toNeighbour
            (
                const tmp<Field<Type>>& fld
            ) const;


        // I/O

            //- Return whether or not the given dictionary contains a
            //  mappedPatchBase specification
            static bool specified(const dictionary& dict);

            //- Write as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPatchBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPatchBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

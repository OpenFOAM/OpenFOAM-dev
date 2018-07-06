/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fvMeshSubset

Description
    Post-processing mesh subset tool.  Given the original mesh and the
    list of selected cells, it creates the mesh consisting only of the
    desired cells, with the mapping list for points, faces, and cells.

    Puts all exposed internal faces into either
    - a user supplied patch
    - a newly created patch "oldInternalFaces"

    - setCellSubset is for small subsets. Uses Maps to minimize memory.
    - setLargeCellSubset is for largish subsets (>10% of mesh).
      Uses labelLists instead.

    - setLargeCellSubset does coupled patch subsetting as well. If it detects
      a face on a coupled patch 'losing' its neighbour it will move the
      face into the oldInternalFaces patch.

    - if a user supplied patch is used it is up to the destination
      patchField to handle exposed internal faces (mapping from face -1).
      If not provided the default is to assign the internalField. All the
      basic patch field types (e.g. fixedValue) will give a warning and
      preferably derived patch field types should be used that know how to
      handle exposed faces (e.g. use uniformFixedValue instead of fixedValue)

SourceFiles
    fvMeshSubset.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshSubset_H
#define fvMeshSubset_H

#include "fvMesh.H"
#include "pointMesh.H"
#include "GeometricField.H"
#include "HashSet.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshSubset Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSubset
{

private:

    // Private data

        //- Mesh to subset from
        const fvMesh& baseMesh_;

        //- Subset mesh pointer
        autoPtr<fvMesh> fvMeshSubsetPtr_;

        //- Point mapping array
        labelList pointMap_;

        //- Face mapping array
        labelList faceMap_;

        //- Cell mapping array
        labelList cellMap_;

        //- Patch mapping array
        labelList patchMap_;

        //- Optional face mapping array with flip encoded
        mutable autoPtr<labelList> faceFlipMapPtr_;


    // Private Member Functions

        //- Check if subset has been performed
        bool checkCellSubset() const;

        //- Mark points in Map
        static void markPoints(const labelList&, Map<label>&);

        //- Mark points (with 0) in labelList
        static void markPoints(const labelList&, labelList&);

        //- Adapt nCellsUsingFace for coupled faces becoming 'uncoupled'.
        void doCoupledPatches
        (
            const bool syncPar,
            Map<label>& facesToSubset,
            labelList& nCellsUsingFace
        ) const;

        //- Subset of subset
        static labelList subset
        (
            const label nElems,
            const labelList& selectedElements,
            const labelList& subsetMap
        );

        //- Create zones for submesh
        void subsetZones();

        //- Helper: extract cells-to-remove from cells-to-keep
        labelList getCellsToRemove
        (
            const labelList& region,
            const label currentRegion
        ) const;

        //- Disallow default bitwise copy construct
        fvMeshSubset(const fvMeshSubset&);

        //- Disallow default bitwise assignment
        void operator=(const fvMeshSubset&);

public:

    // Constructors

        //- Construct given a mesh to subset
        explicit fvMeshSubset(const fvMesh&);


    // Member Functions

        // Edit

            //- Set the subset. Create "oldInternalFaces" patch for exposed
            //  internal faces (patchID==-1) or use supplied patch.
            //  Does not handle coupled patches correctly if only one side
            //  gets deleted.
            void setCellSubset
            (
                const labelHashSet& globalCellMap,
                const label patchID = -1,
                const bool syncPar = true
            );

            //- Set the subset from all cells with region == currentRegion.
            //  Create "oldInternalFaces" patch for exposed
            //  internal faces (patchID==-1) or use supplied patch.
            //  Handles coupled patches by if necessary making coupled patch
            //  face part of patchID (so uncoupled)
            void setLargeCellSubset
            (
                const labelList& region,
                const label currentRegion,
                const label patchID = -1,
                const bool syncCouples = true
            );

            //- setLargeCellSubset but with labelHashSet.
            void setLargeCellSubset
            (
                const labelHashSet& globalCellMap,
                const label patchID = -1,
                const bool syncPar = true
            );


            //- Two step subsetting

                //- Get labels of exposed faces.
                //  These are
                //  - internal faces that become boundary faces
                //  - coupled faces that become uncoupled (since one of the
                //    sides gets deleted)
                labelList getExposedFaces
                (
                    const labelList& region,
                    const label currentRegion,
                    const bool syncCouples = true
                ) const;

                //- For every exposed face (from above getExposedFaces)
                //  used supplied (existing!) patch
                void setLargeCellSubset
                (
                    const labelList& region,
                    const label currentRegion,
                    const labelList& exposedFaces,
                    const labelList& patchIDs,
                    const bool syncCouples = true
                );


        // Access

            //- Original mesh
            const fvMesh& baseMesh() const
            {
                return baseMesh_;
            }

            //- Have subMesh?
            bool hasSubMesh() const;

            //- Return reference to subset mesh
            const fvMesh& subMesh() const;

            fvMesh& subMesh();

            //- Return point map
            const labelList& pointMap() const;

            //- Return face map
            const labelList& faceMap() const;

            //- Return face map with sign to encode flipped faces
            const labelList& faceFlipMap() const;

            //- Return cell map
            const labelList& cellMap() const;

            //- Return patch map
            const labelList& patchMap() const;


        // Field mapping

            //- Map volume field
            template<class Type>
            static tmp<GeometricField<Type, fvPatchField, volMesh>>
            interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                const fvMesh& sMesh,
                const labelList& patchMap,
                const labelList& cellMap,
                const labelList& faceMap
            );

            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh>>
            interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>&
            ) const;

            //- Map surface field. Optionally negates value if flipping
            //  a face (from exposing an internal face)
            template<class Type>
            static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
            interpolate
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&,
                const fvMesh& sMesh,
                const labelList& patchMap,
                const labelList& cellMap,
                const labelList& faceMap,
                const bool negateIfFlipped = true
            );

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
            interpolate
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&,
                const bool negateIfFlipped = true
            ) const;

            //- Map point field
            template<class Type>
            static tmp<GeometricField<Type, pointPatchField, pointMesh>>
            interpolate
            (
                const GeometricField<Type, pointPatchField, pointMesh>&,
                const pointMesh& sMesh,
                const labelList& patchMap,
                const labelList& pointMap
            );

            template<class Type>
            tmp<GeometricField<Type, pointPatchField, pointMesh>>
            interpolate
            (
                const GeometricField<Type, pointPatchField, pointMesh>&
            ) const;

            //- Map dimensioned field
            template<class Type>
            static tmp<DimensionedField<Type, volMesh>>
            interpolate
            (
                const DimensionedField<Type, volMesh>&,
                const fvMesh& sMesh,
                const labelList& cellMap
            );

            template<class Type>
            tmp<DimensionedField<Type, volMesh>>
            interpolate(const DimensionedField<Type, volMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshSubsetInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::fvMeshDistribute

Description
    Sends/receives parts of mesh+fvfields to neighbouring processors.
    Used in load balancing.

    Input is per local cell the processor it should move to. Moves meshes
    and volFields/surfaceFields and returns map which can be used to
    distribute other.

    Notes:
    - does not handle cyclics. Will probably handle separated proc patches.
    - if all cells move off processor also all its processor patches will
      get deleted so comms might be screwed up (since e.g. globalMeshData
      expects procPatches on all)
    - initial mesh has to have procPatches last and all normal patches common
      to all processors and in the same order. This is checked.
    - faces are matched topologically but points on the faces are not. So
      expect problems -on separated patches (cyclics?) -on zero sized processor
      edges.

SourceFiles
    fvMeshDistribute.C
    fvMeshDistributeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshDistribute_H
#define fvMeshDistribute_H

#include "fvMeshSubset.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapAddedPolyMesh;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                      Class fvMeshDistribute Declaration
\*---------------------------------------------------------------------------*/

class fvMeshDistribute
{
    // Private Data

        //- Underlying fvMesh
        fvMesh& mesh_;


    // Private Member Functions

        static void inplaceRenumberWithFlip
        (
            const labelUList& oldToNew,
            const bool oldToNewHasFlip,
            const bool lstHasFlip,
            labelUList& lst
        );

        //- Find indices with value
        static labelList select
        (
            const bool selectEqual,
            const labelList& values,
            const label value
        );

        wordList fieldNames(const word& typeName, label& nFields) const;

        //- Check all procs have same names and in exactly same order
        //  and return the number of fields
        static label checkEqualWordList(const word& typeName, const wordList&);

        //- Merge wordlists over all processors
        static wordList mergeWordList(const wordList&);


        // Patch handling

            //- Find internal patch to put exposed internal faces into
            //  Used when fields are present and require mapping
            label findInternalPatch() const;

            //- Find non-empty patch to put exposed internal faces into
            //  Used when no fields require mapping
            label findNonEmptyPatch() const;

            //- Save boundary fields
            template<class Type, class Mesh>
            void saveBoundaryFields
            (
                PtrList<FieldField<surfaceMesh::PatchField, Type>>&
            ) const;

            //- Map boundary fields
            template<class Type, class Mesh>
            void mapBoundaryFields
            (
                const polyTopoChangeMap& map,
                const PtrList<FieldField<surfaceMesh::PatchField, Type>>&
            );

            //- Set value of exposed patch faces
            template<class Type>
            void initMapExposedFaces(PtrList<Field<Type>>& iflds) const;

            //- Set value of exposed patch faces
            template<class Type>
            void mapExposedFaces
            (
                const polyTopoChangeMap& map,
                const PtrList<Field<Type>>& oldFlds
            );

            //- Correct coupled patch fields
            template<class GeoField>
            void correctCoupledPatchFields();

            //- Delete all processor patches. Move any processor faces into
            //  destinationPatch
            autoPtr<polyTopoChangeMap> deleteProcPatches
            (
                const label destinationPatch
            );

            //- Repatch the mesh. This is only necessary for the proc
            //  boundary faces. newPatchID is over all boundary faces: -1 or
            //  new patchID. constructFaceMap is being adapted for the
            //  possible new face position (since proc faces get automatically
            //  matched)
            autoPtr<polyTopoChangeMap> repatch
            (
                const labelList& newPatchID,
                labelListList& constructFaceMap
            );

            //- Merge any local points that were remotely coupled.
            //  constructPointMap is adapted for the new point labels.
            //
            //  Detect shared points. Need processor patches to be present.
            //  Background: when adding bits of mesh one can get points which
            //  share the same position but are only detectable to be
            //  topologically the same point when doing parallel analysis. This
            //  routine will merge those points.
            autoPtr<polyTopoChangeMap> mergeSharedPoints
            (
                const labelList& pointToGlobalMaster,
                labelListList& constructPointMap
            );


        // Coupling information

            //- Construct the local environment of all boundary faces.
            void getCouplingData
            (
                const labelList& distribution,
                labelList& sourceFace,
                labelList& sourceProc,
                labelList& sourcePatch,
                labelList& sourceNbrPatch,
                labelList& sourceNewProc,
                labelList& sourcePointMaster
            ) const;

            // Subset the neighbourCell/neighbourProc fields
            static void subsetCouplingData
            (
                const fvMesh& mesh,
                const labelList& pointMap,
                const labelList& faceMap,
                const labelList& cellMap,

                const labelList& oldDistribution,
                const labelList& oldFaceOwner,
                const labelList& oldFaceNeighbour,
                const label oldInternalFaces,

                const labelList& sourceFace,
                const labelList& sourceProc,
                const labelList& sourcePatch,
                const labelList& sourceNbrPatch,
                const labelList& sourceNewProc,
                const labelList& sourcePointMaster,

                labelList& subFace,
                labelList& subProc,
                labelList& subPatch,
                labelList& subNbrPatch,
                labelList& subNewProc,
                labelList& subPointMaster
            );

            //- Find cells on mesh whose faceID/procID match the neighbour
            //  cell/proc of domainMesh. Store the matching face.
            static void findCouples
            (
                const primitiveMesh&,
                const labelList& sourceFace,
                const labelList& sourceProc,
                const labelList& sourcePatch,

                const label domain,
                const primitiveMesh& domainMesh,
                const labelList& domainFace,
                const labelList& domainProc,
                const labelList& domainPatch,

                labelList& masterCoupledFaces,
                labelList& slaveCoupledFaces
            );

            //- Map data on boundary faces to new mesh (resulting from adding
            //  two meshes)
            static labelList mapBoundaryData
            (
                const primitiveMesh& mesh,      // mesh after adding
                const mapAddedPolyMesh& map,
                const labelList& boundaryData0, // mesh before adding
                const label nInternalFaces1,
                const labelList& boundaryData1  // added mesh
            );

            //- Map point data to new mesh (resulting from adding two meshes)
            static labelList mapPointData
            (
                const primitiveMesh& mesh,      // mesh after adding
                const mapAddedPolyMesh& map,
                const labelList& boundaryData0, // on mesh before adding
                const labelList& boundaryData1  // on added mesh
            );

        // Other

            //- Remove cells. Add all exposed faces to patch oldInternalPatchi
            autoPtr<polyTopoChangeMap> doRemoveCells
            (
                const labelList& cellsToRemove,
                const label oldInternalPatchi
            );

            //- Add processor patches. Changes mesh and returns per neighbour
            //  proc the processor patchID.
            void addProcPatches
            (
                const labelList& nbrProc,         // Processor that neighbour
                                                  // is now on
                const labelList& referPatchID,    // Original patch ID (or -1)
                const labelList& referNbrPatchID, // Original neighbour patch
                                                  // ID (or -1)
                List<Map<label>>& procPatchID
            );

            //- Add non-conformal processor cyclic patches. Changes the mesh.
            void addNccProcPatches();

            //- Get boundary faces to be repatched. Is -1 or new patchID
            static labelList getBoundaryPatch
            (
                const labelList& neighbourNewProc,  // new processor per b. face
                const labelList& referPatchID,      // -1 or original patch
                const List<Map<label>>& procPatchID // patchID
            );

            //- Send mesh and coupling data.
            static void sendMesh
            (
                const label domain,
                const fvMesh& mesh,
                const wordList& pointZoneNames,
                const wordList& facesZoneNames,
                const wordList& cellZoneNames,
                const labelList& sourceFace,
                const labelList& sourceProc,
                const labelList& sourcePatch,
                const labelList& sourceNbrPatch,
                const labelList& sourceNewProc,
                const labelList& sourcePointMaster,
                Ostream& toDomain
            );

            //- Send subset of fields
            template<class GeoField>
            static void sendFields
            (
                const label domain,
                const wordList& fieldNames,
                const fvMeshSubset&,
                Ostream& toNbr
            );

            //- Receive mesh. Opposite of sendMesh
            static autoPtr<fvMesh> receiveMesh
            (
                const label domain,
                const wordList& pointZoneNames,
                const wordList& facesZoneNames,
                const wordList& cellZoneNames,
                const Time& runTime,
                labelList& domainSourceFace,
                labelList& domainSourceProc,
                labelList& domainSourcePatch,
                labelList& domainSourceNbrPatch,
                labelList& domainSourceNewProc,
                labelList& domainSourcePointMaster,
                Istream& fromNbr
            );

            //- Receive fields. Opposite of sendFields
            template<class GeoField>
            static void receiveFields
            (
                const label domain,
                const wordList& fieldNames,
                typename GeoField::Mesh&,
                PtrList<GeoField>&,
                const dictionary& fieldDicts
            );


public:

    ClassName("fvMeshDistribute");


    // Constructors

        //- Construct from mesh
        fvMeshDistribute(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        fvMeshDistribute(const fvMeshDistribute&) = delete;


    // Member Functions

        //- Helper function: count cells per processor in wanted distribution
        static labelList countCells(const labelList&);

        //- Send cells to neighbours according to distribution
        //  (for every cell the new proc)
        autoPtr<polyDistributionMap> distribute(const labelList& dist);

        // Debugging

            //- Print some info on coupling data
            static void printCoupleInfo
            (
                const primitiveMesh&,
                const labelList&,
                const labelList&,
                const labelList&,
                const labelList&
            );

            //- Print some field info
            template<class GeoField>
            static void printFieldInfo(const fvMesh&);

            //- Print some info on mesh.
            static void printMeshInfo(const fvMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshDistribute&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshDistributeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

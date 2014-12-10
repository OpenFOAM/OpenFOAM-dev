/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::globalPoints

Description
    Calculates points shared by more than two processor patches or cyclic
    patches.

    Is used in globalMeshData. (this info is needed for point/edge
    communication where processor swaps are not enough to exchange data)

    Works purely topological and using local communication only.
    Needs:
      - domain to be one single domain (i.e. all faces can be reached through
        face-cell walk).
      - patch face ordering to be ok
      - f[0] ordering on patch faces to be ok.

    Works by constructing equivalence lists for all the points on processor
    patches. These list are in globalIndexAndTransform numbering
    E.g.
    \verbatim
          ((7 93)(4 731)(3 114))
    \endverbatim

    means point 93 on proc7 is connected to point 731 on proc4 and 114 on proc3.
    It then assigns the lowest numbered processor to be the local 'master' and
    constructs a mapDistribute to send all data to this master.

    Algorithm:
        - get meshPoints of all my points on processor patches and initialize
          equivalence lists to this.
     loop
        - send to all neighbours in relative form:
            - patchFace
            - index in face
        - receive and convert into meshPoints. Add to to my equivalence lists.
        - mark meshPoints for which information changed.
        - send data for these meshPoints again
     endloop until nothing changes

    At this point one will have complete point-point connectivity for all
    points on processor patches. Now (optionally) remove point
    equivalences of size 2. These are just normal points shared
    between two neighbouring procPatches.

    Note: the data held is either mesh point labels (construct from mesh only)
    or patch point labels (construct from mesh and patch).

SourceFiles
    globalPoints.C

\*---------------------------------------------------------------------------*/

#ifndef globalPoints_H
#define globalPoints_H

#include "DynamicList.H"
#include "indirectPrimitivePatch.H"
#include "globalIndex.H"
#include "globalIndexAndTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyBoundaryMesh;
class cyclicPolyPatch;
class polyPatch;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                           Class globalPoints Declaration
\*---------------------------------------------------------------------------*/

class globalPoints
{
    // Private data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Global numbering of untransformed points
        globalIndex globalIndices_;

        //- Global numbering of transformed points
        const globalIndexAndTransform globalTransforms_;

        //- Sum of points on processor patches (unfiltered, point on 2 patches
        //  counts as 2)
        const label nPatchPoints_;

        //- All points on boundaries and their corresponding connected points
        //  on other processors.
        DynamicList<labelPairList> procPoints_;

        //- Map from mesh (or patch) point to index in procPoints
        Map<label> meshToProcPoint_;


        // Calculated mapDistribute addressing

            //- Non-transformed connected points per point (in mapDistribute
            //  indices)
            labelListList pointPoints_;

            //- Transformed points per point (in mapDistribute indices)
            labelListList transformedPointPoints_;

            //- Corresponding map
            autoPtr<mapDistribute> map_;



    // Private Member Functions

        //- Count all points on processorPatches. Is all points for which
        //  information is collected.
        static label countPatchPoints(const polyBoundaryMesh&);

        //- Find index of same processor+index
        label findSamePoint
        (
            const labelPairList& allInfo,
            const labelPair& info
        ) const;

        labelPairList addSendTransform
        (
            const label patchI,
            const labelPairList& info
        ) const;

        //- Add information about patchPointI in relative indices to send
        //  buffers (patchFaces, indexInFace etc.)
        void addToSend
        (
            const polyPatch&,
            const label patchPointI,
            const labelPairList&,
            DynamicList<label>& patchFaces,
            DynamicList<label>& indexInFace,
            DynamicList<labelPairList>& allInfo
        ) const;

        //- Merge info from neighbour into my data
        bool mergeInfo
        (
            const labelPairList& nbrInfo,
            const label localPointI,
            labelPairList& myInfo
        ) const;

        //- From mesh point to 'local point'. Is the mesh point itself
        //  if meshToPatchPoint is empty.
        static label meshToLocalPoint
        (
            const Map<label>& meshToPatchPoint,
            const label meshPointI
        );

        //- Opposite of meshToLocalPoint.
        static label localToMeshPoint
        (
            const labelList& patchToMeshPoint,
            const label localPointI
        );

        //- Store (and merge) info for meshPointI
        bool storeInitialInfo
        (
            const labelPairList& nbrInfo,
            const label localPointI
        );

        //- Store (and merge) info for meshPointI
        bool mergeInfo
        (
            const labelPairList& nbrInfo,
            const label localPointI
        );

        //- Debug printing
        void printProcPoint
        (
            const labelList& patchToMeshPoint,
            const labelPair& pointInfo
        ) const;

        void printProcPoints
        (
            const labelList& patchToMeshPoint,
            const labelPairList& pointInfo
        ) const;

        //- Initialize procPoints_ to my patch points. allPoints = true:
        //  seed with all patch points, = false: only boundaryPoints().
        void initOwnPoints
        (
            const Map<label>& meshToPatchPoint,
            const bool allPoints,
            labelHashSet& changedPoints
        );

        //- Send subset of procPoints to neighbours
        void sendPatchPoints
        (
            const bool mergeSeparated,
            const Map<label>&,
            PstreamBuffers&,
            const labelHashSet&
        ) const;

        //- Receive neighbour points and merge into my procPoints.
        void receivePatchPoints
        (
            const bool mergeSeparated,
            const Map<label>&,
            const labelList&,
            PstreamBuffers&,
            labelHashSet&
        );

        //- Remove entries of size 2 where meshPoint is in provided Map.
        //  Used to remove normal face-face connected points.
        void remove(const labelList& patchToMeshPoint, const Map<label>&);

        //- Return mesh points of other side in same order as my meshPoints.
        static labelList reverseMeshPoints(const cyclicPolyPatch&);

        //- Do all calculations.
        void calculateSharedPoints
        (
            const Map<label>&,
            const labelList&,
            const bool keepAllPoints,
            const bool mergeSeparated
        );

        //- Disallow default bitwise copy construct
        globalPoints(const globalPoints&);

        //- Disallow default bitwise assignment
        void operator=(const globalPoints&);


public:

        //- Declare name of the class and its debug switch
        ClassName("globalPoints");


    // Constructors

        //- Construct from mesh.
        //  keepAllPoints = false : filter out points that are on two
        //  neighbouring coupled patches only (so can be swapped)
        //  mergeSeparated:
        //      true  : merge coupled points across separated patches.
        //      false : do not merge across coupled separated patches.
        globalPoints
        (
            const polyMesh& mesh,
            const bool keepAllPoints,
            const bool mergeSeparated
        );

        //- Construct from mesh and patch of coupled faces. Difference with
        //  construct from mesh only is that this stores the meshToProcPoint,
        //  procPoints as patch local point labels instead of mesh point labels.
        globalPoints
        (
            const polyMesh& mesh,
            const indirectPrimitivePatch& coupledPatch,
            const bool keepAllPoints,
            const bool mergeSeparated
        );


    // Member Functions

        // Access

            //- Global numbering of untransformed (mesh or patch) points
            const globalIndex& globalIndices() const
            {
                return globalIndices_;
            }

            //- Global numbering of transformed (mesh or patch) points
            const globalIndexAndTransform& globalTransforms() const
            {
                return globalTransforms_;
            }

            //- Non-transformed connected points per point (in mapDistribute
            //  indices)
            const labelListList& pointPoints() const
            {
                return pointPoints_;
            }

            //- Non-transformed connected points per point (in mapDistribute
            //  indices)
            labelListList& pointPoints()
            {
                return pointPoints_;
            }

            //- Transformed points per point (in mapDistribute indices)
            const labelListList& transformedPointPoints() const
            {
                return transformedPointPoints_;
            }

            //- Transformed points per point (in mapDistribute indices)
            labelListList& transformedPointPoints()
            {
                return transformedPointPoints_;
            }

            //- Corresponding map
            const mapDistribute& map() const
            {
                return map_();
            }

            //- Corresponding map
            mapDistribute& map()
            {
                return map_();
            }

            //- From (mesh or patch) point to index in procPoints
            const Map<label>& meshToProcPoint() const
            {
                return meshToProcPoint_;
            }

            //- procPoints is per point the connected points (in
            //  globalTransformAndIndex point numbers)
            const DynamicList<labelPairList>& procPoints() const
            {
                return procPoints_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

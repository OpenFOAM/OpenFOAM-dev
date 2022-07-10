/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::polyMesh

Description
    Mesh consisting of general polyhedral cells.

SourceFiles
    polyMesh.C
    polyMeshInitMesh.C
    polyMeshClear.C
    polyMeshFromShapeMesh.C
    polyMeshIO.C
    polyMeshUpdate.C
    polyMeshCheck.C

\*---------------------------------------------------------------------------*/

#ifndef polyMesh_H
#define polyMesh_H

#include "objectRegistry.H"
#include "primitiveMesh.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"
#include "cellShapeList.H"
#include "pointIOField.H"
#include "faceIOList.H"
#include "labelIOList.H"
#include "polyBoundaryMesh.H"
#include "boundBox.H"
#include "meshPointZones.H"
#include "meshFaceZones.H"
#include "meshCellZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class globalMeshData;
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;
class polyMeshTetDecomposition;
class treeDataCell;
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                          Class polyMesh Declaration
\*---------------------------------------------------------------------------*/

class polyMesh
:
    public objectRegistry,
    public primitiveMesh
{

public:

    // Public data types

        //- Enumeration defining the state of the mesh after a read update.
        //  Used for post-processing applications, where the mesh
        //  needs to update based on the files written in time
        //  directories
        enum readUpdateState
        {
            UNCHANGED,
            POINTS_MOVED,
            TOPO_CHANGE,
            TOPO_PATCH_CHANGE
        };

        //- Enumeration defining the decomposition of the cell for
        //  inside/outside test
        enum cellDecomposition
        {
            FACE_PLANES,      //- Faces considered as planes

            FACE_CENTRE_TRIS, //- Faces decomposed into triangles
                              //  using face-centre

            FACE_DIAG_TRIS,   //- Faces decomposed into triangles diagonally

            CELL_TETS         //- Cell decomposed into tets
        };


private:

    // Permanent data

        // Primitive mesh data

            //- Points
            pointIOField points_;

            //- Faces
            faceCompactIOList faces_;

            //- Face owner
            labelIOList owner_;

            //- Face neighbour
            labelIOList neighbour_;

            //- Have the primitives been cleared
            bool clearedPrimitives_;


            //- Boundary mesh
            mutable polyBoundaryMesh boundary_;

            //- Mesh bounding-box.
            //  Created from points on construction, updated when the mesh moves
            boundBox bounds_;

            //- Communicator used for parallel communication
            label comm_;

            //- Vector of non-constrained directions in mesh
            //  defined according to the presence of empty and wedge patches
            mutable Vector<label> geometricD_;

            //- Vector of valid directions in mesh
            //  defined according to the presence of empty patches
            mutable Vector<label> solutionD_;

            //- Base point for face decomposition into tets
            mutable autoPtr<labelIOList> tetBasePtIsPtr_;

            //- Search tree to allow spatial cell searching
            mutable autoPtr<indexedOctree<treeDataCell>> cellTreePtr_;


        // Zoning information

            //- Point zones
            meshPointZones pointZones_;

            //- Face zones
            meshFaceZones faceZones_;

            //- Cell zones
            meshCellZones cellZones_;


        //- Parallel info
        mutable autoPtr<globalMeshData> globalMeshDataPtr_;


        // Mesh motion related dat

            //- Current time index for mesh motion
            mutable label curMotionTimeIndex_;

            //- Old points (for the last mesh motion)
            mutable autoPtr<pointField> oldPointsPtr_;

            //- Old cell centres (for the last mesh motion)
            mutable autoPtr<pointField> oldCellCentresPtr_;

            //- Whether or not to store the old cell centres
            mutable bool storeOldCellCentres_;


    // Private Member Functions

        //- Initialise the polyMesh from the primitive data
        void initMesh();

        //- Initialise the polyMesh from the given set of cells
        void initMesh(cellList& c);

        //- Calculate the valid directions in the mesh from the boundaries
        void calcDirections() const;

        //- Calculate the cell shapes from the primitive
        //  polyhedral information
        void calcCellShapes() const;

        //- Read and return the tetBasePtIs
        autoPtr<labelIOList> readTetBasePtIs() const;

        //- Set the write option of the points
        void setPointsWrite(const IOobject::writeOption wo);

        //- Set the write option of the topology
        void setTopologyWrite(const IOobject::writeOption wo);


        // Helper functions for constructor from cell shapes

            labelListList cellShapePointCells(const cellShapeList&) const;

            labelList facePatchFaceCells
            (
                const faceList& patchFaces,
                const labelListList& pointCells,
                const faceListList& cellsFaceShapes,
                const label patchID
            ) const;

            void setTopology
            (
                const cellShapeList& cellsAsShapes,
                const faceListList& boundaryFaces,
                const wordList& boundaryPatchNames,
                labelList& patchSizes,
                labelList& patchStarts,
                label& defaultPatchStart,
                label& nFaces,
                cellList& cells
            );


        // Geometry checks

            //- Check non-orthogonality
            bool checkFaceOrthogonality
            (
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                const bool detailedReport,
                labelHashSet* setPtr
            ) const;

            //- Check face skewness
            bool checkFaceSkewness
            (
                const pointField& points,
                const vectorField& fCtrs,
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                const bool detailedReport,
                labelHashSet* setPtr
            ) const;

            bool checkEdgeAlignment
            (
                const pointField& p,
                const bool report,
                const Vector<label>& directions,
                labelHashSet* setPtr
            ) const;

            bool checkCellDeterminant
            (
                const vectorField& faceAreas,
                const bool report,
                labelHashSet* setPtr,
                const Vector<label>& meshD
            ) const;

            bool checkFaceWeight
            (
                const vectorField& fCtrs,
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                const scalar minWeight,
                labelHashSet* setPtr
            ) const;

            bool checkVolRatio
            (
                const scalarField& cellVols,
                const bool report,
                const scalar minRatio,
                labelHashSet* setPtr
            ) const;


protected:

    // Protected Member Data

        //- Member data pending transfer to fvMesh

            //- Is the mesh moving
            bool moving_;

            //- Is the mesh topology changing
            bool topoChanged_;


public:

    // Public Typedefs

        typedef polyMesh Mesh;
        typedef polyBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("polyMesh");

    //- Return the default region name
    static word defaultRegion;

    //- Return the mesh sub-directory name (usually "polyMesh")
    static word meshSubDir;


    // Constructors

        //- Construct from IOobject
        explicit polyMesh(const IOobject& io);

        //- Move construct from IOobject or from components.
        //  Boundary is added using addPatches() member function
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& owner,
            labelList&& neighbour,
            const bool syncPar = true
        );

        //- Move construct without boundary with cells rather than
        //  owner/neighbour.
        //  Boundary is added using addPatches() member function
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );

        //- Move construct from cell shapes
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const wordList& boundaryPatchTypes,
            const word& defaultBoundaryPatchName,
            const word& defaultBoundaryPatchType,
            const wordList& boundaryPatchPhysicalTypes,
            const bool syncPar = true
        );

        //- Move construct from cell shapes with patch information in dictionary
        //  format.
        polyMesh
        (
            const IOobject& io,
            pointField&& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const PtrList<dictionary>& boundaryDicts,
            const word& defaultBoundaryPatchName,
            const word& defaultBoundaryPatchType,
            const bool syncPar = true
        );

        //- Disallow default bitwise copy construction
        polyMesh(const polyMesh&);


    //- Destructor
    virtual ~polyMesh();


    // Member Functions

        // Database

            //- Override the objectRegistry dbDir for a single-region case
            virtual const fileName& dbDir() const;

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return the current instance directory for points
            //  Used in the construction of gemometric mesh data dependent
            //  on points
            const fileName& pointsInstance() const;

            //- Return the current instance directory for faces
            const fileName& facesInstance() const;

            //- Set the instance for mesh files
            void setInstance(const fileName&);

            //- Set the instance for the points files
            void setPointsInstance(const fileName&);


        // Access

            // Primitive mesh data

                //- Return raw points
                virtual const pointField& points() const;

                //- Return raw faces
                virtual const faceList& faces() const;

                //- Return face owner
                virtual const labelList& faceOwner() const;

                //- Return face neighbour
                virtual const labelList& faceNeighbour() const;

                //- Return old points for mesh motion
                virtual const pointField& oldPoints() const;

                //- Return old cell centres for mesh motion
                virtual const pointField& oldCellCentres() const;


            //- Return true if io is up-to-date with points
            bool upToDatePoints(const regIOobject& io) const;

            //- Set io to be up-to-date with points
            void setUpToDatePoints(regIOobject& io) const;

            //- Return boundary mesh
            const polyBoundaryMesh& boundaryMesh() const
            {
                return boundary_;
            }

            //- Return mesh bounding box
            const boundBox& bounds() const
            {
                return bounds_;
            }

            //- Return the vector of geometric directions in mesh.
            //  Defined according to the presence of empty and wedge patches.
            //  1 indicates unconstrained direction and -1 a constrained
            //  direction.
            const Vector<label>& geometricD() const;

            //- Return the number of valid geometric dimensions in the mesh
            label nGeometricD() const;

            //- Return the vector of solved-for directions in mesh.
            //  Differs from geometricD in that it includes for wedge cases
            //  the circumferential direction in case of swirl.
            //  1 indicates valid direction and -1 an invalid direction.
            const Vector<label>& solutionD() const;

            //- Return the number of valid solved-for dimensions in the mesh
            label nSolutionD() const;

            //- Return the tetBasePtIs
            const labelIOList& tetBasePtIs() const;

            //- Return the cell search tree
            const indexedOctree<treeDataCell>& cellTree() const;

            //- Return point zones
            const meshPointZones& pointZones() const
            {
                return pointZones_;
            }

            //- Return face zones
            const meshFaceZones& faceZones() const
            {
                return faceZones_;
            }

            //- Return cell zones
            const meshCellZones& cellZones() const
            {
                return cellZones_;
            }

            //- Return parallel info
            const globalMeshData& globalData() const;

            //- Return communicator used for parallel communication
            label comm() const;

            //- Return communicator used for parallel communication
            label& comm();

            //- Return the object registry
            const objectRegistry& thisDb() const
            {
                return *this;
            }


        // Mesh motion

            //- Is mesh moving
            bool moving() const
            {
                return moving_;
            }

            //- Is mesh topology changing
            bool topoChanged() const
            {
                return topoChanged_;
            }

            //- Is mesh changing (topology changing and/or moving)
            bool changing() const
            {
                return moving() || topoChanged();
            }

            //- Reset the points
            //  without storing old points or returning swept volumes
            virtual void setPoints(const pointField&);

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> movePoints(const pointField&);

            //- Reset motion
            void resetMotion() const;


        // Topological change

            //- Return non-const access to the pointZones
            meshPointZones& pointZones()
            {
                return pointZones_;
            }

            //- Return non-const access to the faceZones
            meshFaceZones& faceZones()
            {
                return faceZones_;
            }

            //- Return non-const access to the cellZones
            meshCellZones& cellZones()
            {
                return cellZones_;
            }

            //- Add boundary patches
            void addPatches
            (
                const List<polyPatch*>&,
                const bool validBoundary = true
            );

            //- Add mesh zones
            void addZones
            (
                const List<pointZone*>& pz,
                const List<faceZone*>& fz,
                const List<cellZone*>& cz
            );

            //- Add/insert single patch. If validBoundary the new situation
            //  is consistent across processors.
            virtual void addPatch
            (
                const label insertPatchi,
                const polyPatch& patch,
                const dictionary& patchFieldDict,
                const word& defaultPatchFieldType,
                const bool validBoundary
            );

            //- Reorder and trim existing patches. If validBoundary the new
            //  situation is consistent across processors
            virtual void reorderPatches
            (
                const labelUList& newToOld,
                const bool validBoundary
            );

            //- Update the mesh based on the mesh files saved in
            //  time directories
            virtual readUpdateState readUpdate();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap& map);

            //- Remove boundary patches
            void removeBoundary();

            //- Reset mesh primitive data. Assumes all patch info correct
            //  (so does e.g. parallel communication). If not use
            //  validBoundary=false
            void resetPrimitives
            (
                pointField&& points,
                faceList&& faces,
                labelList&& owner,
                labelList&& neighbour,
                const labelList& patchSizes,
                const labelList& patchStarts,
                const bool validBoundary = true
            );

            //- Reset mesh
            //  For run-time mesh replacement and mesh to mesh mapping
            void reset(const polyMesh&);


        // Storage management

            //- Clear geometry
            void clearGeom();

            //- Clear addressing
            void clearAddressing(const bool isMeshUpdate = false);

            //- Clear all geometry and addressing unnecessary for CFD
            void clearOut();

            //- Clear primitive data (points, faces and cells)
            void clearPrimitives();

            //- Clear tet base points
            void clearTetBasePtIs();

            //- Clear cell tree data
            void clearCellTree();

            //- Remove all files from mesh instance
            void removeFiles(const fileName& instanceDir) const;

            //- Remove all files from mesh instance()
            void removeFiles() const;


        // Geometric checks. Selectively override primitiveMesh functionality.

            //- Check non-orthogonality
            virtual bool checkFaceOrthogonality
            (
                const bool report = false,
                labelHashSet* setPtr = nullptr
            ) const;

            //- Check face skewness
            virtual bool checkFaceSkewness
            (
                const bool report = false,
                labelHashSet* setPtr = nullptr
            ) const;

            //- Check edge alignment for 1D/2D cases
            virtual bool checkEdgeAlignment
            (
                const bool report,
                const Vector<label>& directions,
                labelHashSet* setPtr
            ) const;

            virtual bool checkCellDeterminant
            (
                const bool report,
                labelHashSet* setPtr
            ) const;

            //- Check mesh motion for correctness given motion points
            virtual bool checkMeshMotion
            (
                const pointField& newPoints,
                const bool report = false,
                const bool detailedReport = false
            ) const;

            //- Check for face weights
            virtual bool checkFaceWeight
            (
                const bool report,
                const scalar minWeight = 0.05,
                labelHashSet* setPtr = nullptr
            ) const;

            //- Check for neighbouring cell volumes
            virtual bool checkVolRatio
            (
                const bool report,
                const scalar minRatio = 0.01,
                labelHashSet* setPtr = nullptr
            ) const;


        // Position search functions

            //- Find the cell, tetFacei and tetPti for point p
            void findCellFacePt
            (
                const point& p,
                label& celli,
                label& tetFacei,
                label& tetPti
            ) const;

            //- Find the tetFacei and tetPti for point p in celli.
            //  tetFacei and tetPtI are set to -1 if not found
            void findTetFacePt
            (
                const label celli,
                const point& p,
                label& tetFacei,
                label& tetPti
            ) const;

            //- Test if point p is in the celli
            bool pointInCell
            (
                const point& p,
                label celli,
                const cellDecomposition = CELL_TETS
            ) const;

            //- Find cell enclosing this location and return index
            //  If not found -1 is returned
            label findCell
            (
                const point& p,
                const cellDecomposition = CELL_TETS
            ) const;


        // Write

            //- Write the underlying polyMesh
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write = true
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

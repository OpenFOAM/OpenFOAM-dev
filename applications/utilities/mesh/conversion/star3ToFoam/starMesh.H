/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::starMesh

Description
    A messy mesh class which supports the possibility of creating a shapeMesh
    for regular Star meshes (no arbitrary interfaces or collapsed SAMM cells).
    If any of these special feateres exist, the mesh is created as polyMesh

SourceFiles
    calcPointCells.C
    createBoundaryFaces.C
    createCoupleMatches.C
    createPolyBoundary.C
    createPolyCells.C
    fixCollapsedEdges.C
    mergeCoupleFacePoints.C
    readBoundary.C
    readCells.C
    readCouples.C
    readPoints.C
    starMesh.C
    writeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef starMesh_H
#define starMesh_H

#include "polyMesh.H"
#include "cellShape.H"
#include "cellList.H"
#include "polyPatchList.H"
#include "coupledFacePair.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class starMesh Declaration
\*---------------------------------------------------------------------------*/

class starMesh
{
    // Private Data

        //- Name of the case
        fileName casePrefix_;

        //- Database
        const Time& runTime_;

        //- Points supporting the mesh
        pointField points_;

        //- Cell shapes
        cellShapeList cellShapes_;

        //- Boundary faces
        faceListList boundary_;

        //- Boundary patch types
        wordList patchTypes_;

        //- Default boundary patch name
        word defaultFacesName_;

        //- Default boundary patch types
        word defaultFacesType_;

        //- Boundary patch names
        wordList patchNames_;

        //- Boundary patch physical types
        wordList patchPhysicalTypes_;

        //- Point labels (STAR point numbering is not necessarily contiguous)
        labelList starPointLabelLookup_;

        //- STAR point number for a given vertex
        labelList starPointID_;

        //- STAR Cell number for a given cell
        labelList starCellID_;

        //- Cell labels (STAR cell numbering is not necessarily contiguous)
        labelList starCellLabelLookup_;

        //- STAR Cell permutation label
        labelList starCellPermutation_;

        //- List of faces for every cell
        faceListList cellFaces_;

        //- Cell ID for every boundary face. Used in two-tier boundary
        //  reconstruction
        labelListList boundaryCellIDs_;

        //- Cell face ID for every boundary face. Used in two-tier boundary
        //  reconstruction
        labelListList boundaryCellFaceIDs_;

        //- Global face list for polyMesh
        faceList meshFaces_;

        //- Cells as polyhedra for polyMesh
        cellList cellPolys_;

        //- Number of internal faces for polyMesh
        label nInternalFaces_;

        //- Polyhedral mesh boundary patch start indices
        labelList polyBoundaryPatchStartIndices_;

        //- Point-cell addressing. Used for topological analysis
        // Warning. This point cell addressing list potentially contains
        // duplicate cell entries. Use additional checking
        mutable labelListList* pointCellsPtr_;

        //- List of face couples
        PtrList<coupledFacePair> couples_;

        //- Can the mesh be treated as a shapeMesh?
        bool isShapeMesh_;

    // Private static data members

        //- Error on unity small tolerance
        static const scalar smallMergeTol_;

        //- Couple match relative point merge tolerance
        static const scalar cpMergePointTol_;

        //- Pointers to cell models
        static const cellModel* unknownPtr_;
        static const cellModel* tetPtr_;
        static const cellModel* pyrPtr_;
        static const cellModel* tetWedgePtr_;
        static const cellModel* prismPtr_;
        static const cellModel* wedgePtr_;
        static const cellModel* hexPtr_;

        static const cellModel* sammTrim1Ptr_;
        static const cellModel* sammTrim2Ptr_;
        static const cellModel* sammTrim3Ptr_;
        static const cellModel* sammTrim4Ptr_;
        static const cellModel* sammTrim5Ptr_;
        static const cellModel* sammTrim8Ptr_;

        //- Regular addressing data
        static const label regularAddressingTable[6][8];

        //- SAMM addressing data
        static const label sammAddressingTable[9][12];

        static const label sammFacePermutationTable[24][8];

        static const label shapeFaceLookup[19][9];


    // Private Member Functions

        //- Read fixed format vertex label
        static label readVtxLabel(IFstream&);

        //- Read fixed format vertex coordinate component
        static scalar readVtxCmpt(IFstream&);

        //- Read to nl
        static void readToNl(IFstream&);

        //- Read the points file
        void readPoints(const scalar scaleFactor);

        //- Read the cells file
        void readCells();

        void addRegularCell
        (
            const labelList& labels,
            const label nCreatedCells
        );

        void addSAMMcell
        (
            const labelList& labels,
            const label nCreatedCells
        );

        //- Read the boundary file
        void readBoundary();

        //- Check and correct collapsed edges on faces
        // Note. If a collapsed edge is found, the mesh is no longer shapeMesh
        void fixCollapsedEdges();

        //- Read couples
        void readCouples();

        //- Create couple matches
        void createCoupleMatches();

        //- Calculate pointCells
        void calcPointCells() const;

        const labelListList& pointCells() const;

        //- Mark boundary faces from the quads
        void markBoundaryFaces();

        //- Collect boundary faces from the quads
        void collectBoundaryFaces();

        //- Specialist version of face comparison to deal with
        // PROSTAR boundary format idiosyncracies
        bool starEqualFace
        (
            const face& boundaryFace,
            const face& cellFace
        ) const;

        //- Merge couple face points
        void mergeCoupleFacePoints();

        //- Point indexing structure used by sort to keep track of
        //  the point labels
        struct pointIndex
        {
            label label_;
            double index_;
        };

        //- Purge cell shapes
        void purgeCellShapes();

        //- Make polyhedral cells and global faces if the mesh is polyhedral
        void createPolyCells();

        //- Make polyhedral boundary from shape boundary
        // (adds more faces to the face list)
        void createPolyBoundary();

        //- Make polyhedral mesh data (packing)
        void createPolyMeshData();

        //- Add polyhedral boundary
        List<polyPatch*> polyBoundaryPatches(const polyMesh&);

        //- Clear extra storage before creation of the mesh to remove
        //  a memory peak
        void clearExtraStorage();


public:

    // Constructors

        //- Construct from case name
        starMesh
        (
            const fileName& prefix,
            const Time& rt,
            const scalar scaleFactor
        );

        //- Disallow default bitwise copy construction
        starMesh(const starMesh&) = delete;


    //- Destructor
    ~starMesh();


    // Member Functions

        //- Write mesh
        void writeMesh();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const starMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

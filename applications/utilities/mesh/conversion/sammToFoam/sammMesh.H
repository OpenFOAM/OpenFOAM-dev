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
    Foam::sammMesh

Description
    A messy mesh class which supports the possibility of creating a shapeMesh
    for regular Samm meshes (no arbitrary interfaces or collapsed SAMM cells).
    If any of these special feateres exist, the mesh is created as polyMesh

SourceFiles
    sammMesh.C

\*---------------------------------------------------------------------------*/

#ifndef sammMesh_H
#define sammMesh_H

#include "polyMesh.H"
#include "cellShape.H"
#include "cellList.H"
#include "polyPatchList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sammMesh Declaration
\*---------------------------------------------------------------------------*/

class sammMesh
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

        //- Point labels (SAMM point numbering is not necessarily contiguous)
        labelList starPointLabelLookup_;

        //- Point labels (SAMM point numbering is not necessarily contiguous)
        labelList starCellLabelLookup_;

        //- List of faces for every cell
        faceListList cellFaces_;

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

        //- Can the mesh be treated as a shapeMesh?
        bool isShapeMesh_;

    // Private static data members

        //- Pointers to cell models
        static const cellModel* unknownPtr_;
        static const cellModel* hexPtr_;
        static const cellModel* wedgePtr_;
        static const cellModel* prismPtr_;
        static const cellModel* pyrPtr_;
        static const cellModel* tetPtr_;
        static const cellModel* tetWedgePtr_;

        static const cellModel* sammTrim1Ptr_;
        static const cellModel* sammTrim2Ptr_;
        static const cellModel* sammTrim3Ptr_;
        static const cellModel* sammTrim4Ptr_;
        static const cellModel* sammTrim5Ptr_;
        static const cellModel* sammTrim8Ptr_;

        static const label shapeFaceLookup[19][9];


        //- SAMM addressing data
        static List<const cellModel*> sammShapeLookup;
        static List<const label*> sammAddressingTable;

    // Private Member Functions

        //- Fill SAMM lookup tables
        void fillSammCellShapeTable();
        void fillSammAddressingTable();


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
            const label typeFlag,
            const labelList& globalLabels,
            const label nCreatedCells
        );


        //- Read the boundary file
        void readBoundary();


        //- Check and correct collapsed edges on faces
        // Note. If a collapsed edge is found, the mesh is no longer shapeMesh
        void fixCollapsedEdges();

        //- Read couples
        void readCouples();

        //- Calculate pointCells
        void calcPointCells() const;

        const labelListList& pointCells() const;

        //- Create boundary faces from the quads
        void createBoundaryFaces();

        //- Specialist version of face comparison to deal with
        // PROSTAR boundary format idiosyncracies
        bool sammEqualFace
        (
            const face& boundaryFace,
            const face& cellFace
        ) const;

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
        List<polyPatch* > polyBoundaryPatches(const polyMesh&);


public:

    // Constructors

        //- Construct from case name
        sammMesh
        (
            const fileName& prefix,
            const Time& rt,
            const scalar scaleFactor
        );

        //- Disallow default bitwise copy construction
        sammMesh(const sammMesh&) = delete;


    //- Destructor
    ~sammMesh();


    // Member Functions

        //- Write mesh
        void writeMesh();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sammMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

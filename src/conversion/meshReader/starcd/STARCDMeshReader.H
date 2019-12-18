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
    Foam::meshReaders::STARCD

Description
    Read pro-STAR vrt/cel/bnd files.
    The protected data in meshReader are filled.

    Starting with pro-STAR version 4, the files have become easier to read.
    - vertices are space-delimited.
    - the cell format is logical.
    - trimmed and degenerate cells are saved as polyhedral.
    - the boundaries corresponds to cells and their faces.

SourceFiles
    STARCDMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDMeshReader_H
#define STARCDMeshReader_H

#include "meshReader.H"
#include "boundaryRegion.H"
#include "cellShape.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace meshReaders
{

/*---------------------------------------------------------------------------*\
                    Class meshReaders::STARCD Declaration
\*---------------------------------------------------------------------------*/

class STARCD
:
    public meshReader
{

protected:

    // Protected Data

        static const char* const defaultBoundaryName;
        static const char* const defaultSolidBoundaryName;

        //- Face addressing from pro-STAR faces -> OpenFOAM faces
        static const int starToFoamFaceAddr[4][6];

        //- Cell shapes
        cellShapeList cellShapes_;

        //- Point labels (imported Point numbering not necessarily contiguous)
        labelList mapToFoamPointId_;

        //- Cell labels (imported Cell numbering not necessarily contiguous)
        labelList mapToFoamCellId_;

        //- Boundary region data
        boundaryRegion boundaryRegion_;


    // Protected Member Functions

        //- Read the mesh from the file(s)
        virtual bool readGeometry(const scalar scaleFactor = 1.0);

        //- Read points from file
        void readPoints(const fileName&, const scalar scaleFactor);

        //- Read cell connectivities from file
        virtual void readCells(const fileName&);

        //- Remove unused points
        void cullPoints();

        //- Read boundary (cell/face) definitions
        void readBoundary(const fileName&);

        //- Read auxiliary data from constant/{boundaryRegion,cellTable}
        void readAux(const objectRegistry&);

        //- Read and discard to newline
        static void readToNewline(IFstream&);

        //- Read header
        static bool readHeader(IFstream&, word fileSignature);


protected:

        enum cellType
        {
            starcdFluidType   = 1,
            starcdSolidType   = 2,
            starcdBaffleType  = 3,
            starcdShellType   = 4,
            starcdLineType    = 5,
            starcdPointType   = 6
        };

        enum shapeType
        {
            starcdPoint = 1,
            starcdLine  = 2,
            starcdShell = 3,
            starcdHex   = 11,
            starcdPrism = 12,
            starcdTet   = 13,
            starcdPyr   = 14,
            starcdPoly  = 255
        };


public:

    // Static data

        //- Keep solids (default false)
        static bool keepSolids;


    // Constructors

        //- Construct from case name
        STARCD
        (
            const fileName& prefix,
            const objectRegistry&,
            const scalar scaleFactor = 1.0
        );

        //- Disallow default bitwise copy construction
        STARCD(const STARCD&) = delete;


    //- Destructor
    virtual ~STARCD();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const STARCD&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshReaders
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

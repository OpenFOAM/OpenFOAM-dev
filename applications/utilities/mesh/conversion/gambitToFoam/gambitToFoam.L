/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    gambitToFoam

Description
    Converts a GAMBIT mesh to OpenFOAM format.

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------- *\
   ------ local definitions
 \* ------------------------------------------------------------------------- */

#include "scalarList.H"
#include "IStringStream.H"

// For EOF only
#include <cstdio>

using namespace Foam;

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "preservePatchTypes.H"
#include "cellModeller.H"
#include "cellShape.H"
#include "SLList.H"
#include "SLPtrList.H"

label nPoints = 0;
label nCells = 0;
label nCellStreams = 0;
label nPatches = 0;

label nCoordDirections = 0;
label nVectorComponents = 0;


pointField points(0);
labelList pointMap(0);

PtrList<labelList> cellLabels(0);
labelList cellMap(0);
labelList cellTypes(0);
labelList cellStreamIDs(0);

wordList patchNames(0);
labelListList patchCells(0);
labelListList patchCellFaces(0);
label nValuesForPatchFaces = 0;

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
#if YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+
spaceNl               ({space}|\n)*

alpha                 [_[:alpha:]]
digit                 [[:digit:]]

dotColonDash          [.:-]

label                 [0-9]{digit}*
zeroLabel             {digit}*

word                  ({alpha}|{digit}|{dotColonDash})*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              ((({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))|0)

x                     {floatNum}
y                     {floatNum}
z                     {floatNum}

labelListElement      {space}{zeroLabel}
scalarListElement     {space}{floatNum}
labelList             ({labelListElement}+{space})
scalarList            ({scalarListElement}+{space})

starStar              ("**")
text                  ({space}({word}*{space})*\n)

dateDDMMYYYY          ({digit}{digit}"/"{digit}{digit}"/"{digit}{digit}{digit}{digit})
dateDDMonYYYY         ((({digit}{digit}{space})|({digit}{space})){alpha}*{space}{digit}{digit}{digit}{digit})
time                  ({digit}{digit}":"{digit}{digit}":"{digit}{digit})

versionNumber         ({digit}|".")*

controlInfo           ^{space}"CONTROL INFO"{space}{versionNumber}
applicationData       ^{space}"APPLICATION DATA"{space}{versionNumber}
nodalCoords           ^{space}"NODAL COORDINATES"{space}{versionNumber}
cellsAndElements      ^{space}"ELEMENTS/CELLS"{space}{versionNumber}
cellStreams           ^{space}"ELEMENT GROUP"{space}{versionNumber}
boundaryPatch         ^{space}"BOUNDARY CONDITIONS"{space}{versionNumber}
faceConnectivity      ^{space}"FACE CONNECTIVITY"{space}{versionNumber}

endOfSection          ^{space}"ENDOFSECTION"{space}

program               {space}"PROGRAM:"{space}
version               {space}"VERSION:"{space}

group                 {space}"GROUP:"{space}
elements              {space}"ELEMENTS:"{space}
material              {space}"MATERIAL:"{space}
nFlags                {space}"NFLAGS:"{space}
mtype                 {space}"MTYPE:"{space}


 /* ------------------------------------------------------------------------- *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------- */

%option stack

%x controlInfo
%x readControlHeader
%x readTitle
%x readProgramID
%x readVersionID
%x applicationData
%x nodalCoords
%x cellsAndElements
%x cellStreams
%x cellStreamTitle
%x cellStreamFlags
%x cellStreamLabels
%x readCellStreamGroupID
%x readCellStreamNElements
%x readCellStreamMaterial
%x readCellStreamNFlags
%x boundaryPatch
%x boundaryPatchParams
%x boundaryPatchFaces
%x faceConnectivity

%x globalMeshData
%x cellContLine

%%

%{
    label nCellContinuationLines = 0;
    label curNumberOfNodes = 0;
    label curNumberOfCells = 0;
    label curGroupID = 0;
    label nFlagsForStream = 0;
    label curBoundaryPatch = 0;
    label curPatchFace = 0;
%}


 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

 /*                      ------ Reading control header ------                 */

{controlInfo} {
        BEGIN(readControlHeader);
    }


<readControlHeader>{starStar}{space}{text} {
        BEGIN(readTitle);
    }


<readTitle>{text} {

        Info  << "Title: " << YYText();
        BEGIN(controlInfo);
    }


<controlInfo>{spaceNl}{program} {
        BEGIN(readProgramID);
    }


<readProgramID>{space}{word} {

        Info<< "Written by " << YYText() << " ";
        BEGIN(controlInfo);
    }


<controlInfo>{spaceNl}{version} {
        BEGIN(readVersionID);
    }


<readVersionID>{space}{versionNumber} {

        Info<< " version " << YYText() << endl;
        BEGIN(controlInfo);
    }


<controlInfo>{space}{dateDDMonYYYY}{space}{time} {

        Info<< "File written on " << YYText() << endl;
    }


<controlInfo>{space}{dateDDMMYYYY}{space}{time} {

        Info<< "File written on " << YYText() << endl;
    }


<controlInfo>{spaceNl}"NUMNP"{space}"NELEM"{space}"NGRPS"{space}"NBSETS"{space}("NDFCD"|"NDCFD"){space}"NDFVL"{space}\n {
        BEGIN(globalMeshData);
    }


<globalMeshData>{spaceNl}{label}{space}{label}{space}{label}{space}{label}{space}{label}{space}{label}{space}\n {

        IStringStream nodeStream(YYText());

        nPoints = readLabel(nodeStream);
        nCells = readLabel(nodeStream);
        nCellStreams = readLabel(nodeStream);
        nPatches = readLabel(nodeStream);
        nCoordDirections = readLabel(nodeStream);
        nVectorComponents = readLabel(nodeStream);

        // reset list sizes - now known!
        points.setSize(nPoints);
        pointMap.setSize(nPoints);

        cellLabels.setSize(nCells);
        cellMap.setSize(nCells);
        cellTypes.setSize(nCells);
        cellStreamIDs.setSize(nCells);

        patchNames.setSize(nPatches);
        patchCells.setSize(nPatches);
        patchCellFaces.setSize(nPatches);

        Info<< "    number of points:  " <<  nPoints << endl
            << "    number of cells:   " << nCells << endl
            << "    number of patches: " << nPatches << endl;

        BEGIN(controlInfo);
    }


 /*                 ------  Reading nodal coordinates ------                 */

{nodalCoords}{spaceNl} {

        curNumberOfNodes = 0;
        Info<< "Reading nodal coordinates" << endl;
        BEGIN(nodalCoords);
    }


<nodalCoords>{spaceNl}{label}{space}{x}{space}{y}{space}{z}{space}\n {

        IStringStream nodeStream(YYText());

        label nodeI(readLabel(nodeStream));

        // Note: coordinates must be read one at the time.
        scalar x = readScalar(nodeStream);
        scalar y = readScalar(nodeStream);
        scalar z = readScalar(nodeStream);

        // add mapping and scalced node to the list
        pointMap[curNumberOfNodes] = nodeI;
        points[curNumberOfNodes] = point(x, y, z);
        curNumberOfNodes++;
    }


 /*                   ------ Reading cells and elements ------               */

{cellsAndElements}{spaceNl} {

        curNumberOfCells = 0;
        Info<< "Reading cells" << endl;
        BEGIN(cellsAndElements);
    }

<cellsAndElements>{spaceNl}{label}{space}{label}{space}{label}{labelList} {

        IStringStream elementStream(YYText());

        label cellI(readLabel(elementStream));
        label cellType(readLabel(elementStream));
        label nVertices(readLabel(elementStream));

        // reset number of continuation lines
        nCellContinuationLines = 0;

        cellMap[curNumberOfCells] = cellI;
        cellTypes[curNumberOfCells] = cellType;
        cellLabels.set(curNumberOfCells, new labelList(nVertices));

        labelList& curLabels = cellLabels[curNumberOfCells];

        // Find out how many labels are expected. If less or equal to
        // seven, read them all and finish with it. If there is more,
        // set read of the next line
        label labelsToRead = min(8, nVertices);
        label labelI = 0;
        for (; labelI < labelsToRead; labelI++)
        {
            if (elementStream.eof()) break;

            // Check token to avoid trailing space.
            token curLabelTok(elementStream);
            if (curLabelTok.isLabel())
            {
                curLabels[labelI] = curLabelTok.labelToken();
            }
            else
            {
                break;
            }
        }

        if (labelI < nVertices)
        {
            BEGIN(cellContLine);
        }
        else
        {
            curNumberOfCells++;
        }
    }


<cellContLine>{spaceNl}{labelList} {

        IStringStream elementStream(YYText());

        nCellContinuationLines++;

        labelList& curLabels = cellLabels[curNumberOfCells];

        label labelsToRead = min
        (
            (nCellContinuationLines + 1)*7,
            curLabels.size()
        );

        for
        (
            label labelI = nCellContinuationLines*7;
            labelI < labelsToRead;
            labelI++
        )
        {
            curLabels[labelI] = readLabel(elementStream);
        }

        // if read is finished, go back to reading cells
        if (curLabels.size() < (nCellContinuationLines + 1)*7)
        {
            curNumberOfCells++;

            BEGIN(cellsAndElements);
        }
    }


 /*             ------ Reading element group information ------               */

{cellStreams}{spaceNl} {
        Info<< "Reading cell streams" << endl;
        BEGIN(cellStreams);
    }


<cellStreams>{spaceNl}{group} {
        BEGIN(readCellStreamGroupID);
    }


<readCellStreamGroupID>{space}{label} {
        IStringStream groupStream(YYText());

        if (curGroupID > 0)
        {
            FatalErrorIn("gambitToFoam::main")
                << "<readCellStreamGroupID>{space}{label} : "
                << "trying to reset group ID while active"
                << abort(FatalError);
        }
        else
        {
            curGroupID = readLabel(groupStream);
        }

        BEGIN(cellStreams);
    }


<cellStreams>{spaceNl}{elements} {
        BEGIN(readCellStreamNElements);
    }


<readCellStreamNElements>{space}{label} {

        IStringStream nElementsStream(YYText());

        readLabel(nElementsStream);

        BEGIN(cellStreams);
    }

<cellStreams>{spaceNl}{material} {
        BEGIN(readCellStreamMaterial);
    }


<readCellStreamMaterial>{space}{label} {

        IStringStream materialIDstream(YYText());

        readLabel(materialIDstream);

        BEGIN(cellStreams);
    }


<cellStreams>{spaceNl}{nFlags} {
        BEGIN(readCellStreamNFlags);
    }


<readCellStreamNFlags>{space}{label} {

        IStringStream nFlagsStream(YYText());

        nFlagsForStream = readLabel(nFlagsStream);

        BEGIN(cellStreamTitle);
    }


<cellStreamTitle>{spaceNl}{word}{spaceNl} {

        word streamName(Foam::string::validate<word>(YYText()));

        BEGIN(cellStreamFlags);
    }


<cellStreamFlags>{labelList} {
        Info<< "Reading cell stream labels" << endl;
        BEGIN(cellStreamLabels);
    }


<cellStreamLabels>{spaceNl}{labelList} {

        IStringStream nFlagsStream(YYText());

        label cellLabel;
        while (nFlagsStream.read(cellLabel))
        {
            cellStreamIDs[cellLabel - 1] = curGroupID;
        }

        // reset current group ID and a number of flags
        curGroupID = 0;
        nFlagsForStream = 0;
    }


 /*            ------  Reading end of section and others ------               */

<cellStreamLabels>{endOfSection}\n {

        Info<< "Finished reading cell stream labels" << endl;

        // reset current group ID and a number of flags
        curGroupID = 0;
        nFlagsForStream = 0;

        BEGIN(INITIAL);
    }


{boundaryPatch}{spaceNl} {
        curPatchFace = 0;
        Info<< "Reading patches" << endl;
        BEGIN(boundaryPatchParams);
    }


<boundaryPatchParams>{spaceNl}{word}{labelList} {

        IStringStream patchParamsStream(YYText());

        patchParamsStream.read(patchNames[curBoundaryPatch]);

        readLabel(patchParamsStream);
        label nEntry(readLabel(patchParamsStream));
        nValuesForPatchFaces = readLabel(patchParamsStream);

        patchCells[curBoundaryPatch].setSize(nEntry);
        patchCellFaces[curBoundaryPatch].setSize(nEntry);

        Info<< "patch " << curBoundaryPatch
            << ": name: " << patchNames[curBoundaryPatch]
            << endl;

        BEGIN(boundaryPatchFaces);
    }


<boundaryPatchFaces>{spaceNl}{label}{space}{label}{space}{label}({scalarList}|{space}\n) {

        // Face-based boundary condition
        IStringStream patchFacesStream(YYText());

        patchCells[curBoundaryPatch][curPatchFace] =
            readLabel(patchFacesStream);
        readLabel(patchFacesStream);
        patchCellFaces[curBoundaryPatch][curPatchFace] =
            readLabel(patchFacesStream);

        // patch face values currently discarded
        if (nValuesForPatchFaces > 0)
        {
            scalarList patchFaceValues(nValuesForPatchFaces);

            forAll(patchFaceValues, fI)
            {
                patchFaceValues[fI] = readScalar(patchFacesStream);
            }
        }

        curPatchFace++;
    }


<boundaryPatchFaces>{spaceNl}{label}({scalarList}|\n) {

        // Vertex-based boundary condition
        FatalErrorIn("gambitToFoam::main")
            << "<boundaryPatchFaces>{spaceNl}{label}{scalarList} : "
            << "boundary condition specified on vertices not supported"
            << abort(FatalError);
    }


<boundaryPatchFaces>{endOfSection}\n {

        curBoundaryPatch++;

        BEGIN(INITIAL);
    }


 /*             ------ Reading end of section and others ------               */

<controlInfo,nodalCoords,cellsAndElements>{endOfSection}\n {
        BEGIN(INITIAL);
    }


 /* ------ Ignore remaining space and \n s.  Any other characters are errors. */

.|\n {}


 /*  ------ On EOF return to previous file, if none exists terminate. ------  */

<<EOF>> {
            yyterminate();
    }
%%


#include "fileName.H"
#include <fstream>
using std::ifstream;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("GAMBIT file");
    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1"
    );

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    const scalar scaleFactor = args.optionLookupOrDefault("scale", 1.0);

#   include "createTime.H"

    const fileName gambitFile = args[1];
    ifstream gambitStream(gambitFile.c_str());

    if (!gambitStream)
    {
        FatalErrorIn("gambitToFoam::main")
            << args.executable()
            << ": file " << gambitFile << " not found"
            << abort(FatalError);
    }

    yyFlexLexer lexer(&gambitStream);
    while (lexer.yylex() != 0)
    {}

    Info<< "Finished lexing" << endl;

    // make a point mapping array
    label maxPointIndex = 0;

    forAll(pointMap, pointI)
    {
        if (pointMap[pointI] > maxPointIndex)
        {
            maxPointIndex = pointMap[pointI];
        }
    }


    labelList pointLookup(maxPointIndex + 1, -1);

    forAll(pointMap, pointI)
    {
        pointLookup[pointMap[pointI] ] = pointI;
    }

    // make a cell mapping array
    label maxCellIndex = 0;

    forAll(cellMap, cellI)
    {
        if (cellMap[cellI] > maxCellIndex)
        {
            maxCellIndex = cellMap[cellI];
        }
    }

    labelList cellLookup(maxCellIndex + 1);

    forAll(cellMap, cellI)
    {
        cellLookup[cellMap[cellI] ] = cellI;
    }

    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    labelList labelsHex(8);
    labelList labelsPrism(6);
    labelList labelsPyramid(5);
    labelList labelsTet(4);

    cellShapeList cells(cellLabels.size());

    forAll(cellTypes, cellI)
    {
        const labelList& curCellLabels = cellLabels[cellI];

        // Tetrahedron
        if (cellTypes[cellI] == 6)
        {
            labelsTet[0] = pointLookup[curCellLabels[0] ];
            labelsTet[1] = pointLookup[curCellLabels[2] ];
            labelsTet[2] = pointLookup[curCellLabels[3] ];
            labelsTet[3] = pointLookup[curCellLabels[1] ];

            cells[cellI] = cellShape(tet, labelsTet);
        }

        // Square-based pyramid
        else if (cellTypes[cellI] == 7)
        {
            labelsPyramid[0] = pointLookup[curCellLabels[0] ];
            labelsPyramid[1] = pointLookup[curCellLabels[1] ];
            labelsPyramid[2] = pointLookup[curCellLabels[3] ];
            labelsPyramid[3] = pointLookup[curCellLabels[2] ];
            labelsPyramid[4] = pointLookup[curCellLabels[4] ];

            cells[cellI] = cellShape(pyr, labelsPyramid);
        }

        // Triangular prism
        else if (cellTypes[cellI] == 5)
        {
            labelsPrism[0] = pointLookup[curCellLabels[0] ];
            labelsPrism[1] = pointLookup[curCellLabels[1] ];
            labelsPrism[2] = pointLookup[curCellLabels[2] ];
            labelsPrism[3] = pointLookup[curCellLabels[3] ];
            labelsPrism[4] = pointLookup[curCellLabels[4] ];
            labelsPrism[5] = pointLookup[curCellLabels[5] ];

            cells[cellI] = cellShape(prism, labelsPrism);
        }

        // Hex
        else if (cellTypes[cellI] == 4)
        {
            labelsHex[0] = pointLookup[curCellLabels[0] ];
            labelsHex[1] = pointLookup[curCellLabels[1] ];
            labelsHex[2] = pointLookup[curCellLabels[3] ];
            labelsHex[3] = pointLookup[curCellLabels[2] ];
            labelsHex[4] = pointLookup[curCellLabels[4] ];
            labelsHex[5] = pointLookup[curCellLabels[5] ];
            labelsHex[6] = pointLookup[curCellLabels[7] ];
            labelsHex[7] = pointLookup[curCellLabels[6] ];

            cells[cellI] = cellShape(hex, labelsHex);
        }
    }

    // give foam model face number given a fluent model face number
    label faceIndex[8][6] =
    {
        {-1, -1, -1, -1, -1, -1},  // 0
        {-1, -1, -1, -1, -1, -1},  // 1
        {-1, -1, -1, -1, -1, -1},  // 2
        { 2,  1,  3,  0,  4,  5},  // Hex (3)
        {-1, -1, -1, -1, -1, -1},  // 4
        { 4,  3,  2,  0,  1, -1},  // Triangular prism (5)
        { 0,  4,  3,  2,  1, -1},  // Pyramid (6)
        { 2,  1,  0,  3, -1, -1}   // Tet (7)
    };

    faceListList boundary(patchCells.size());

    forAll(patchCells, patchI)
    {
        labelList& curCells = patchCells[patchI];
        labelList& curFaces = patchCellFaces[patchI];

        faceList& patchFaces = boundary[patchI];
        patchFaces.setSize(curCells.size());

        forAll(curCells, faceI)
        {
            patchFaces[faceI] =
                cells[cellLookup[curCells[faceI] ] ].faces()
                [
                    faceIndex
                    [
                        // this picks a cell type
                        cells[cellLookup[curCells[faceI] ] ]
                            .model().index()
                    ]
                    [curFaces[faceI] - 1] // this gives a fluent face - 1
                ];

        }
    }

    Info<< "gambitToFoam: " << endl
        << "Gambit file format does not provide information about the type of "
        << "the patch (eg. wall, symmetry plane, cyclic etc)." << endl
        << "All the patches have been created "
        << "as type patch. Please reset after mesh conversion as necessary."
        << endl;

    // Scale points
    points *= scaleFactor;

    PtrList<dictionary> patchDicts(boundary.size());
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );

    // Add information to dictionary
    forAll(patchNames, patchI)
    {
        if (!patchDicts.set(patchI))
        {
            patchDicts.set(patchI, new dictionary());
        }
        // Add but not overwrite
        patchDicts[patchI].add("type", polyPatch::typeName, false);
    }

    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove(points),
        cells,
        boundary,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< "Writing polyMesh" << endl;
    pShapeMesh.write();

    Info<< "\nEnd\n" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of gambitToFoam.L
 \* ------------------------------------------------------------------------- */

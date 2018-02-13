/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

Application
    fluentMeshToFoam

Description
    Converts a Fluent mesh to OpenFOAM format
    including multiple region and region boundary handling.

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------- *\
   ------ local definitions
 \* ------------------------------------------------------------------------- */

#include "argList.H"
#include "Time.H"
#include "IStringStream.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cellShape.H"
#include "faceSet.H"
#include "cellSet.H"
#include "meshTools.H"
#include "OFstream.H"

#include "readHexLabel.H"
#include "cellShapeRecognition.H"
#include "repatchPolyTopoChanger.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

const scalar convertToMeters = 1.0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


label dimensionOfGrid = 0;

label nPoints = 0;
label nFaces = 0;
label nCells = 0;

pointField points(0);
SLList<label> pointGroupZoneID;
SLList<label> pointGroupStartIndex;
SLList<label> pointGroupEndIndex;


faceList faces(0);
labelList owner(0);
labelList neighbour(0);
SLList<label> faceGroupZoneID;
SLList<label> faceGroupStartIndex;
SLList<label> faceGroupEndIndex;

labelList fluentCellModelID(0);
SLList<label> cellGroupZoneID;
SLList<label> cellGroupStartIndex;
SLList<label> cellGroupEndIndex;
SLList<label> cellGroupType;

// number of zones adjusted at run-time if necessary
label maxZoneID = 100;
label zoneIDBuffer = 10;

wordList patchTypeIDs(maxZoneID);
wordList patchNameIDs(maxZoneID);

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
#if YY_FLEX_MINOR_VERSION < 6 && YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}

%}

one_space                  [ \t\f]
space                      {one_space}*
some_space                 {one_space}+
cspace                     ","{space}
spaceNl                    ({space}|\n|\r)*

alpha                      [_[:alpha:]]
digit                      [[:digit:]]
decDigit                   [[:digit:]]
octalDigit                 [0-7]
hexDigit                   [[:xdigit:]]

lbrac                      "("
rbrac                      ")"
quote                      \"
dash                       "-"
dotColonDash               [.:-]

schemeSpecialInitial       [!$%&*/:<=>?~_^#.]
schemeSpecialSubsequent    [.+-]
schemeSymbol               (({some_space}|{alpha}|{quote}|{schemeSpecialInitial})({alpha}|{quote}|{digit}|{schemeSpecialInitial}|{schemeSpecialSubsequent})*)


identifier                 {alpha}({alpha}|{digit})*
integer                    {decDigit}+
label                      [1-9]{decDigit}*
hexLabel                   {hexDigit}+
zeroLabel                  {digit}*

word                       ({alpha}|{digit}|{dotColonDash})*

exponent_part              [eE][-+]?{digit}+
fractional_constant        [-+]?(({digit}*"."{digit}+)|({digit}+".")|({digit}))

double                     ((({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))|0)

x                          {double}
y                          {double}
z                          {double}
scalar                     {double}
labelListElement           {space}{zeroLabel}
hexLabelListElement        {space}{hexLabel}
scalarListElement          {space}{double}
schemeSymbolListElement    {space}{schemeSymbol}
labelList                  ({labelListElement}+{space})
hexLabelList               ({hexLabelListElement}+{space})
scalarList                 ({scalarListElement}+{space})
schemeSymbolList           ({schemeSymbolListElement}+{space})

starStar                   ("**")
text                       ({space}({word}*{space})*)

dateDDMMYYYY               ({digit}{digit}"/"{digit}{digit}"/"{digit}{digit}{digit}{digit})
dateDDMonYYYY              ((({digit}{digit}{space})|({digit}{space})){alpha}*{space}{digit}{digit}{digit}{digit})
time                       ({digit}{digit}":"{digit}{digit}":"{digit}{digit})

versionNumber              ({digit}|".")*

comment                    {spaceNl}"(0"{space}
header                     {spaceNl}"(1"{space}
dimension                  {spaceNl}"(2"{space}
point                      {spaceNl}"(10"{space}
fluentFace                 {spaceNl}"(13"{space}
cell                       {spaceNl}"(12"{space}
zoneVariant1               {spaceNl}"(39"{space}
zoneVariant2               {spaceNl}"(45"{space}

unknownPeriodicFace        {spaceNl}"(17"{space}
periodicFace               {spaceNl}"(18"{space}
cellTree                   {spaceNl}"(58"{space}
faceTree                   {spaceNl}"(59"{space}
faceParents                {spaceNl}"(61"{space}

endOfSection               {space}")"{space}



 /* ------------------------------------------------------------------------- *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------- */

%option stack

%x readComment
%x embeddedCommentState
%x readHeader
%x readDimension
%x readPoint
%x readPointHeader
%x readNumberOfPoints
%x readPointGroupData
%x readPointData
%x readPoints2D
%x readPoints3D
%x fluentFace
%x readFaceHeader
%x readNumberOfFaces
%x readFaceGroupData
%x readFaceData
%x readFacesMixed
%x readFacesUniform
%x cell
%x readCellHeader
%x readNumberOfCells
%x readCellGroupData
%x readCellData
%x readCellsMixed
%x readCellsUniform
%x zone
%x readZoneHeader
%x readZoneGroupData
%x readZoneData
%x readZoneBlock

%x periodicFace
%x cellTree
%x faceTree
%x faceParents

%x unknownBlock
%x embeddedUnknownBlock
%%

%{
    // Point data
    label pointGroupNumberOfComponents = 3;
    label pointi = 0; // index used for reading points

    // Face data
    label faceGroupElementType = -1;
    label facei = 0;

    // Cell data
    label cellGroupElementType = -1;
    label celli = 0;
%}


 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

 /*                      ------ Reading control header ------                 */

{comment} {
        yy_push_state(readComment);
    }


<readComment>{quote}{text}{quote} {
    }


<readComment>{spaceNl}{endOfSection} {
        yy_pop_state();
    }

{header} {
        BEGIN(readHeader);
    }

<readHeader>{quote}{text}{quote} {
        Info<< "Reading header: " << YYText() << endl;
    }


{dimension} {
        BEGIN(readDimension);
    }

<readDimension>{space}{label}{space} {
        IStringStream dimOfGridStream(YYText());

        dimensionOfGrid = readLabel(dimOfGridStream);

        Info<< "Dimension of grid: " << dimensionOfGrid << endl;
    }


{point} {
        BEGIN(readPointHeader);
    }

<readPointHeader>{spaceNl}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfPoints);
    }

<readNumberOfPoints>{space}{hexLabel}{space}{labelList} {

        IStringStream numberOfPointsStream(YYText());

        nPoints = readHexLabel(numberOfPointsStream);

        Info<< "Number of points: " << nPoints << endl;

        // set the size of the points list
        points.setSize(nPoints);

        // meaningless type skipped
        readLabel(numberOfPointsStream);

        // this dimension of grid may be checked against global dimension
        if (numberOfPointsStream)
        {
            // check dimension of grid
            readLabel(numberOfPointsStream);
        }
        else
        {
            Info<< endl;
        }
    }

<readPointHeader>{spaceNl}{lbrac} {
        BEGIN(readPointGroupData);
    }

<readPointGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{labelList} {
        IStringStream pointGroupDataStream(YYText());

        // read point zone-ID, start and end-label
        // the indices will be used for checking later.
        pointGroupZoneID.append(readHexLabel(pointGroupDataStream));

        pointGroupStartIndex.append(readHexLabel(pointGroupDataStream));

        pointGroupEndIndex.append(readHexLabel(pointGroupDataStream));

        // point group type skipped
        readHexLabel(pointGroupDataStream);

        // In FOAM, indices start from zero - adjust
        pointi = pointGroupStartIndex.last() - 1;

        // reset number of components to default
        pointGroupNumberOfComponents = 3;

        // read number of components in the vector
        if (pointGroupDataStream)
        {
            pointGroupNumberOfComponents = readLabel(pointGroupDataStream);
        }
    }

<readNumberOfPoints,readPointGroupData>{endOfSection} {
        BEGIN(readPointData);
    }

<readPointData>{spaceNl}{lbrac} {

        Info<< "Reading points" << endl;

        if (pointGroupNumberOfComponents == 2)
        {
            yy_push_state(readPoints2D);
        }
        else
        {
            yy_push_state(readPoints3D);
        }
    }

<readPoints2D>{spaceNl}{scalarList} {

        IStringStream vertexXyzStream(YYText());

        // Note: coordinates must be read one at the time.
        scalar x = readScalar(vertexXyzStream);
        scalar y = readScalar(vertexXyzStream);

        points[pointi] = point(x, y, 0);
        pointi++;
    }

<readPoints3D>{spaceNl}{scalarList} {

        IStringStream vertexXyzStream(YYText());

        // Note: coordinates must be read one at the time.
        scalar x = readScalar(vertexXyzStream);
        scalar y = readScalar(vertexXyzStream);
        scalar z = readScalar(vertexXyzStream);

        points[pointi] = convertToMeters*point(x, y, z);
        pointi++;
    }

<readPoints2D,readPoints3D>{spaceNl}{endOfSection} {

        // check read of points
        if (pointi != pointGroupEndIndex.last())
        {
            Info<< "problem with reading points: "
                << "start index: " << pointGroupStartIndex.last()
                << " end index: " << pointGroupEndIndex.last()
                << " last points read: " << pointi << endl;
        }

        yy_pop_state();
    }


{fluentFace} {
        BEGIN(readFaceHeader);
    }

<readFaceHeader>{spaceNl}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfFaces);
    }

<readNumberOfFaces>{space}{hexLabel}{space}{labelListElement}+ {

        IStringStream numberOfFacesStream(YYText());

        nFaces = readHexLabel(numberOfFacesStream);

        Info<< "number of faces: " << nFaces << endl;

        faces.setSize(nFaces);
        owner.setSize(nFaces);
        neighbour.setSize(nFaces);

        // Meaningless type and element type not read
    }

<readFaceHeader>{spaceNl}{lbrac} {
        BEGIN(readFaceGroupData);
    }

<readFaceGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{hexLabelListElement}+ {

        IStringStream faceGroupDataStream(YYText());

        // read fluentFace zone-ID, start and end-label
        faceGroupZoneID.append(readHexLabel(faceGroupDataStream));

        // the indices will be used for checking later.
        faceGroupStartIndex.append(readHexLabel(faceGroupDataStream));

        faceGroupEndIndex.append(readHexLabel(faceGroupDataStream));

        // face group type
        readHexLabel(faceGroupDataStream);

        faceGroupElementType = readHexLabel(faceGroupDataStream);

        // In FOAM, indices start from zero - adjust
        facei = faceGroupStartIndex.last() - 1;
    }

<readNumberOfFaces,readFaceGroupData>{spaceNl}{endOfSection} {
        BEGIN(readFaceData);
    }

<readFaceData>{spaceNl}{lbrac} {

        if (faceGroupElementType == 0)
        {
            Info<< "Reading mixed faces" << endl;
            yy_push_state(readFacesMixed);
        }
        else
        {
            Info<< "Reading uniform faces" << endl;
            yy_push_state(readFacesUniform);
        }
    }

<readFacesMixed>{spaceNl}{hexLabelList} {

        IStringStream mixedFaceStream(YYText());

        face& curFaceLabels = faces[facei];

        // set size of label list
        curFaceLabels.setSize(readLabel(mixedFaceStream));

        forAll(curFaceLabels, i)
        {
            curFaceLabels[i] = readHexLabel(mixedFaceStream) - 1;
        }

        // read neighbour and owner. Neighbour comes first
        neighbour[facei] = readHexLabel(mixedFaceStream) - 1;
        owner[facei] = readHexLabel(mixedFaceStream) - 1;
        facei++;
    }

<readFacesUniform>{spaceNl}{hexLabelList} {

        IStringStream mixedFaceStream(YYText());

        face& curFaceLabels = faces[facei];

        // set size of label list. This is OK because in Fluent the type
        // for edge is 2, for triangle is 3 and for quad is 4
        curFaceLabels.setSize(faceGroupElementType);

        forAll(curFaceLabels, i)
        {
            curFaceLabels[i] = readHexLabel(mixedFaceStream) - 1;
        }

        // read neighbour and owner. Neighbour comes first
        neighbour[facei] = readHexLabel(mixedFaceStream) - 1;
        owner[facei] = readHexLabel(mixedFaceStream) - 1;
        facei++;
    }

<readFacesMixed,readFacesUniform>{spaceNl}{endOfSection} {

        // check read of fluentFaces
        if (facei != faceGroupEndIndex.last())
        {
            Info<< "problem with reading fluentFaces: "
                << "start index: " << faceGroupStartIndex.last()
                << " end index: " << faceGroupEndIndex.last()
                << " last fluentFaces read: " << facei << endl;
        }

        yy_pop_state();
    }


{cell} {
        BEGIN(readCellHeader);
    }

<readCellHeader>{spaceNl}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfCells);
    }

<readNumberOfCells>{space}{hexLabel}{space}{labelListElement}+ {

        IStringStream numberOfCellsStream(YYText());

        nCells = readHexLabel(numberOfCellsStream);

        Info<< "Number of cells: " << nCells << endl;

        fluentCellModelID.setSize(nCells);

        // Meaningless type and element type not read
    }

<readCellHeader>{spaceNl}{lbrac} {
        BEGIN(readCellGroupData);
    }

<readCellGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel} {
        // Warning. This entry must be above the next one because of the lexing
        // rules. It is introduced to deal with the problem of reading
        // non-standard cell definition from Tgrid, which misses the type label.

        Info<< "Tgrid syntax problem: " << YYText() << endl;
        IStringStream cellGroupDataStream(YYText());

        // read cell zone-ID, start and end-label
        cellGroupZoneID.append(readHexLabel(cellGroupDataStream));

        // the indices will be used for checking later.
        cellGroupStartIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupEndIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupType.append(readHexLabel(cellGroupDataStream));

        Info<< "cellGroupZoneID:" << cellGroupZoneID.last()
            << endl;
        Info<< "cellGroupStartIndex:" << cellGroupStartIndex.last()
            << endl;
        Info<< "cellGroupEndIndex:" << cellGroupEndIndex.last()
            << endl;
        Info<< "cellGroupType:" << cellGroupType.last()
            << endl;


        // Note. Potentially skip cell set if type is zero.
        // This entry does not exist in Tgrid files.
        if (dimensionOfGrid == 2)
        {
            // Tgrid creating triangles
            cellGroupElementType = 1;
        }
        else
        {
            cellGroupElementType = 2;
        }

        // In FOAM, indices start from zero - adjust
        celli = cellGroupStartIndex.last() - 1;

        if (cellGroupElementType != 0)
        {
            label lastIndex = cellGroupEndIndex.last();

            for (; celli < lastIndex; celli++)
            {
                fluentCellModelID[celli] = cellGroupElementType;
            }
        }
    }

<readCellGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel} {
        // Warning. See above

        Info<< "Other readCellGroupData: " << YYText() << endl;


        IStringStream cellGroupDataStream(YYText());

        // read cell zone-ID, start and end-label
        cellGroupZoneID.append(readHexLabel(cellGroupDataStream));

        // the indices will be used for checking later.
        cellGroupStartIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupEndIndex.append(readHexLabel(cellGroupDataStream));

        cellGroupType.append(readHexLabel(cellGroupDataStream));

        // Note. Potentially skip cell set if type is zero.

        cellGroupElementType = readHexLabel(cellGroupDataStream);

        // In FOAM, indices start from zero - adjust
        celli = cellGroupStartIndex.last() - 1;

        if (cellGroupElementType != 0)
        {
            Info<< "Reading uniform cells" << endl;
            label lastIndex = cellGroupEndIndex.last();

            for (; celli < lastIndex; celli++)
            {
                fluentCellModelID[celli] = cellGroupElementType;
            }
        }
    }

<readNumberOfCells,readCellGroupData>{endOfSection} {
        BEGIN(readCellData);
    }

<readCellData>{spaceNl}{lbrac} {
        Info<< "Reading mixed cells" << endl;
        yy_push_state(readCellsMixed);
    }

<readCellsMixed>{spaceNl}{labelList} {

        IStringStream fluentCellModelIDStream(YYText());

        label celliD;
        while (fluentCellModelIDStream.read(celliD))
        {
            fluentCellModelID[celli] = celliD;
            celli++;
        }
    }

<readCellsMixed>{spaceNl}{endOfSection} {

        // check read of cells
        if (celli != cellGroupEndIndex.last())
        {
            Info<< "problem with reading cells: "
                << "start index: " << cellGroupStartIndex.last()
                << " end index: " << cellGroupEndIndex.last()
                << " last cells read: " << celli << endl;
        }

        yy_pop_state();
    }



{zoneVariant1} {
        BEGIN(readZoneHeader);
    }

{zoneVariant2} {
        BEGIN(readZoneHeader);
    }

<readZoneHeader>{spaceNl}{lbrac} {
        BEGIN(readZoneGroupData);
    }

<readZoneGroupData>{space}{label}{space}{word}{space}{word} {

        IStringStream zoneDataStream(YYText());

        // cell zone-ID not in hexadecimal!!! Inconsistency
        label zoneID(readLabel(zoneDataStream));

        if (zoneID > maxZoneID - 1)
        {
            // resize the container
            maxZoneID = zoneID + zoneIDBuffer;

            patchTypeIDs.setSize(maxZoneID);
            patchNameIDs.setSize(maxZoneID);
        }

        zoneDataStream >> patchTypeIDs[zoneID];
        zoneDataStream >> patchNameIDs[zoneID];

        Info<< "Read zone1:" << zoneID
            << " name:" << patchNameIDs[zoneID]
            << " patchTypeID:" << patchTypeIDs[zoneID]
            << endl;
    }

<readZoneGroupData>{space}{label}{space}{word}{space}{word}{space}{label} {
        // Fluent manual inconsistency, version 6.1.22
        IStringStream zoneDataStream(YYText());

        // cell zone-ID not in hexadecimal!!! Inconsistency
        label zoneID(readLabel(zoneDataStream));

        if (zoneID > maxZoneID - 1)
        {
            // resize the container
            maxZoneID = zoneID + zoneIDBuffer;

            patchTypeIDs.setSize(maxZoneID);
            patchNameIDs.setSize(maxZoneID);
        }

        zoneDataStream >> patchTypeIDs[zoneID];
        zoneDataStream >> patchNameIDs[zoneID];

        Info<< "Read zone2:" << zoneID
            << " name:" << patchNameIDs[zoneID]
            << " patchTypeID:" << patchTypeIDs[zoneID]
            << endl;
    }

<readZoneGroupData>{endOfSection} {
        BEGIN(readZoneData);
    }

<readZoneData>{spaceNl}{lbrac} {
        Info<< "Reading zone data" << endl;
        yy_push_state(readZoneBlock);
    }

<readZoneBlock>{spaceNl}{schemeSymbolList} {
    }

<readZoneBlock>{lbrac} {
        Info<< "Found unknown block in zone:" << YYText() << endl;
        yy_push_state(unknownBlock);
    }

<readZoneBlock>{endOfSection} {
        yy_pop_state();
    }



 /*             ------ Reading end of section and others ------               */

<readHeader,readDimension,readPointData,readFaceData,readCellData,readZoneData>{spaceNl}{endOfSection} {
        BEGIN(INITIAL);
    }

 /*    ------ Reading unknown type or non-standard comment ------             */


{lbrac}{label} {
        Info<< "Found unknown block:" << YYText() << endl;
        yy_push_state(unknownBlock);
    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{schemeSymbol} {
    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{lbrac} {
        Info<< "Embedded blocks in comment or unknown:" << YYText() << endl;
        yy_push_state(embeddedUnknownBlock);

    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{endOfSection} {
        Info<< "Found end of section in unknown:" << YYText() << endl;
        yy_pop_state();
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{labelList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{hexLabelList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{scalarList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{schemeSymbolList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{text} {
    }


 /* ------ Ignore remaining space and \n s.  Any other characters are errors. */

<readPoints2D,readPoints3D>.|\n {
    Info<< "ERROR! Do not understand characters: " << YYText() << endl;
    }
.|\n {}


 /*  ------ On EOF return to previous file, if none exists terminate. ------  */

<<EOF>> {
            yyterminate();
    }
%%


#include "fileName.H"
#include <fstream>

// Find label of face.
label findFace(const primitiveMesh& mesh, const face& f)
{
    // Get faces using zeroth vertex of face.
    const labelList& pFaces = mesh.pointFaces()[f[0] ];

    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        if (f == mesh.faces()[facei])
        {
            return facei;
        }
    }

    // Didn't find face. Do what?
    FatalErrorInFunction
        << "Problem : cannot find a single face in the mesh which uses"
        << " vertices " << f << exit(FatalError);

    return -1;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("Fluent mesh file");
    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1"
    );
    argList::addOption
    (
        "2D",
        "thickness",
        "use when converting a 2-D mesh (applied before scale)"
    );
    argList::addBoolOption
    (
        "writeSets",
        "write cell zones and patches as sets"
    );
    argList::addBoolOption
    (
        "writeZones",
        "write cell zones as zones"
    );

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    const scalar scaleFactor = args.optionLookupOrDefault("scale", 1.0);

    const bool writeSets  = args.optionFound("writeSets");
    const bool writeZones = args.optionFound("writeZones");

    #include "createTime.H"

    const fileName fluentFile = args[1];
    std::ifstream fluentStream(fluentFile.c_str());

    if (!fluentStream)
    {
        FatalErrorInFunction
            << args.executable()
            << ": file " << fluentFile << " not found"
            << exit(FatalError);
    }

    yyFlexLexer lexer(&fluentStream);
    while (lexer.yylex() != 0)
    {}

    Info<< "\n\nFINISHED LEXING\n\n\n";
    // Lookup table giving number of vertices given a fluent cell type ID
    // Currently not used.
    // label fluentModelNVertices[7] = {-1, 3, 4, 4, 8, 5, 6};

    // Lookup table giving number of vertices given a fluent cell type ID
    label fluentModelNFaces[7] = {-1, 3, 4, 4, 6, 5, 5};

    // Make a list of cell faces to be filled in for owner and neighbour

    labelListList cellFaces(nCells);

    labelList nFacesInCell(nCells, 0);

    forAll(cellFaces, celli)
    {
        cellFaces[celli].setSize(fluentModelNFaces[fluentCellModelID[celli] ]);
    }

    // fill in owner and neighbour

    forAll(owner, facei)
    {
        if (owner[facei] > -1)
        {
            label curCell = owner[facei];
            cellFaces[curCell][nFacesInCell[curCell] ] = facei;

            nFacesInCell[curCell]++;
        }
    }

    forAll(neighbour, facei)
    {
        if (neighbour[facei] > -1)
        {
            label curCell = neighbour[facei];
            cellFaces[curCell][nFacesInCell[curCell] ] = facei;

            nFacesInCell[curCell]++;
        }
    }

    // Construct shapes from face lists
    cellShapeList cellShapes(nCells);

    Info<< "dimension of grid: " << dimensionOfGrid << endl;
    faceList frontAndBackFaces;

    if (dimensionOfGrid == 2)
    {
        // Extrude 2-D mesh into 3-D, in z-direction

        scalar twoDThickness = 1.0;

        if (!args.optionReadIfPresent("2D", twoDThickness))
        {
            const scalar extrusionFactor = 0.02; //0.01 in each direction
            boundBox box(points);
            twoDThickness = extrusionFactor*box.mag();
        }

        Info<< "Grid is 2-D. Extruding in z-direction by: " << twoDThickness
            << endl;

        const scalar zOffset = twoDThickness / 2;

        pointField oldPoints = points;

        const label pointOffset = oldPoints.size();

        points.setSize(2*pointOffset);

        label nNewPoints = 0;

        // Note: In order for the owner-neighbour rules to be right, the
        // points given by Fluent need to represent the FRONT plane of the
        // geometry. Therefore, the extrusion will be in -z direction
        //
        forAll(oldPoints, pointi)
        {
            points[nNewPoints] = oldPoints[pointi];

            points[nNewPoints].z() = zOffset;

            nNewPoints++;
        }

        forAll(oldPoints, pointi)
        {
            points[nNewPoints] = oldPoints[pointi];

            points[nNewPoints].z() = -zOffset;

            nNewPoints++;
        }

        // 2-D shape recognition
        Info<< "Creating shapes for 2-D cells"<< endl;

        // Set the number of empty faces
        frontAndBackFaces.setSize(2*nCells);

        forAll(fluentCellModelID, celli)
        {
            switch (fluentCellModelID[celli])
            {
                case 1:
                {
                    cellShapes[celli] =
                        cellShape
                        (
                            extrudedTriangleCellShape
                            (
                                celli,
                                cellFaces[celli],
                                faces,
                                owner,
                                neighbour,
                                pointOffset,
                                frontAndBackFaces
                            )
                        );
                }
                break;

                case 3:
                {
                    cellShapes[celli] =
                        cellShape
                        (
                            extrudedQuadCellShape
                            (
                                celli,
                                cellFaces[celli],
                                faces,
                                owner,
                                neighbour,
                                pointOffset,
                                frontAndBackFaces
                            )
                        );
                }
                break;

                default:
                {
                    FatalErrorInFunction
                        << fluentCellModelID[celli]
                        << abort(FatalError);
                }
            }
        }

        // Create new faces
        forAll(faces, facei)
        {

            if (faces[facei].size() != 2)
            {
                FatalErrorInFunction
                    << "fluentMeshToFoam: a 2-D face defined with "
                    << faces[facei].size() << " points." << endl;
            }

            labelList& newFace = faces[facei];

            newFace.setSize(4);

            newFace[2] = newFace[1] + pointOffset;

            newFace[3] = newFace[0] + pointOffset;
        }

        // Create new cells from 2-D shapes
    }
    else
    {
        // 3-D shape recognition
        Info<< "Creating shapes for 3-D cells"<< endl;
        forAll(fluentCellModelID, celli)
        {
            if
            (
                fluentCellModelID[celli] == 2    // tet
             || fluentCellModelID[celli] == 4    // hex
             || fluentCellModelID[celli] == 5    // pyramid
             || fluentCellModelID[celli] == 6    // prism
            )
            {
                cellShapes[celli] =
                    cellShape
                    (
                        create3DCellShape
                        (
                            celli,
                            cellFaces[celli],
                            faces,
                            owner,
                            neighbour,
                            fluentCellModelID[celli]
                        )
                    );
            }
            else
            {
                FatalErrorInFunction
                    << "unrecognised 3-D cell shape: "
                    << fluentCellModelID[celli]
                    << abort(FatalError);
            }
        }
    }

    // boundary faces are oriented such that the owner is zero and the face
    // area vector points into the domain. Turn them round before making patches
    // for Foam compatibility

    forAll(faces, facei)
    {
        if (owner[facei] == -1)
        {
            // reverse face
            labelList oldFace = faces[facei];

            forAllReverse(oldFace, i)
            {
                faces[facei][oldFace.size() - i - 1] =
                    oldFace[i];
            }
        }
    }


    //make patchless mesh before analysing boundaries

    faceListList patches(0);
    wordList patchNames(0);
    wordList patchTypes(0);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    wordList patchPhysicalTypes(0);

    // Scale the points

    points *= scaleFactor;

    Info<< "Building patch-less mesh..." << flush;

    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove(points),
        cellShapes,
        patches,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    // don't write mesh yet, otherwise preservePatchTypes will be broken
    // and zones won't be written
    // checkmesh done after patch addition as well
    Info<< "done." << endl;


    Info<< endl << "Building boundary and internal patches." << endl;
    //adding patches after mesh construction allows topological checks
    //on whether a patch is internal or external, something fluent
    //doesn't seem to mind

    // Make boundary patches

    SLList<label>::iterator faceGroupZoneIDIter = faceGroupZoneID.begin();
    SLList<label>::iterator faceGroupStartIndexIter =
        faceGroupStartIndex.begin();
    SLList<label>::iterator faceGroupEndIndexIter = faceGroupEndIndex.begin();

    // Note. Not all groups of faces will be boundary patches.
    // Take care on construction

    //2D needs extra space for frontAndBack faces
    if (dimensionOfGrid == 2)
    {
        patches.setSize(faceGroupZoneID.size()+1);
        patchNames.setSize(faceGroupZoneID.size()+1);
        patchTypes.setSize(faceGroupZoneID.size()+1);
        patchPhysicalTypes.setSize(faceGroupZoneID.size()+1);
    }
    else
    {
        patches.setSize(faceGroupZoneID.size());
        patchNames.setSize(faceGroupZoneID.size());
        patchTypes.setSize(faceGroupZoneID.size());
        patchPhysicalTypes.setSize(faceGroupZoneID.size());
    }

    label nPatches = 0;

    //colate information for all patches (internal and external)
    for
    (
        ;
        faceGroupZoneIDIter != faceGroupZoneID.end()
     && faceGroupStartIndexIter != faceGroupStartIndex.end()
     && faceGroupEndIndexIter != faceGroupEndIndex.end();
        ++faceGroupZoneIDIter,
        ++faceGroupStartIndexIter,
        ++faceGroupEndIndexIter
    )
    {
        // get face type and name
        const word& curPatchType = patchTypeIDs[faceGroupZoneIDIter()];

        const word& curPatchName = patchNameIDs[faceGroupZoneIDIter()];

        Info<< "Creating patch " << nPatches
            << " for zone: " << faceGroupZoneIDIter()
            << " start: " << faceGroupStartIndexIter()
            << " end: " << faceGroupEndIndexIter()
            << " type: " << curPatchType << " name: " << curPatchName << endl;

        // make patch labels
        label faceLabel = faceGroupStartIndexIter() - 1;

        faceList patchFaces(faceGroupEndIndexIter() - faceLabel);

        forAll(patchFaces, facei)
        {
            if
            (
                faces[faceLabel].size() == 3
             || faces[faceLabel].size() == 4
            )
            {
                patchFaces[facei] = face(faces[faceLabel]);
                faceLabel++;
            }
            else
            {
                FatalErrorInFunction
                    << "unrecognised face shape with "
                    << patchFaces[facei].size() << " vertices"
                    << abort(FatalError);
            }
        }

        //inlets and outlets
        if
        (
            curPatchType == "pressure"
         || curPatchType == "pressure-inlet"
         || curPatchType == "inlet-vent"
         || curPatchType == "intake-fan"
         || curPatchType == "pressure-outlet"
         || curPatchType == "exhaust-fan"
         || curPatchType == "outlet-vent"
         || curPatchType == "pressure-far-field"
         || curPatchType == "velocity-inlet"
         || curPatchType == "mass-flow-inlet"
         || curPatchType == "outflow"
        )
        {
            patches[nPatches] = patchFaces;
            patchTypes[nPatches] = polyPatch::typeName;
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else if (curPatchType == "wall" ) //wall boundaries
        {
            patches[nPatches] = patchFaces;
            patchTypes[nPatches] = wallPolyPatch::typeName;
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else if
        (
            curPatchType == "symmetry"
         || curPatchType == "axis"
        ) //symmetry planes
        {
            patches[nPatches] = patchFaces;
            patchTypes[nPatches] = symmetryPolyPatch::typeName;
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else if
        (
            curPatchType == "interior"
         || curPatchType == "interface"
         || curPatchType == "internal"
         || curPatchType == "solid"
         || curPatchType == "fan"
         || curPatchType == "radiator"
         || curPatchType == "porous-jump"
        ) //interior boundaries - will not be added as patches
        {
            patches[nPatches] = patchFaces;
            patchTypes[nPatches] = "internal";
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else if
        (
            curPatchType == ""
        ) //unnamed face regions default to interior patches
        {
            Info<< "Patch " << faceGroupZoneIDIter()
                << ": Faces are defined but "
                << "not created as a zone." << endl
                << "Null specification is only valid for internal faces."
                << endl;

            patches[nPatches] = patchFaces;
            patchTypes[nPatches] = "internal";
            patchNames[nPatches] = curPatchName;

            nPatches++;
        }
        else //unknown face regions are not handled
        {
            FatalErrorInFunction
                << "fluent patch type " << curPatchType << " not recognised."
                << abort(FatalError);
        }
    }

    //add front and back boundaries for 2D meshes
    if (dimensionOfGrid == 2)
    {
        Info<< "Creating patch for front and back planes" << endl << endl;

        patches[nPatches] = frontAndBackFaces;
        patchTypes[nPatches] = emptyPolyPatch::typeName;
        patchNames[nPatches] = "frontAndBackPlanes";

        nPatches++;
    }

    //Now have all patch information,
    //check whether each patch is internal or external
    //and add boundaries to mesh
    //also write face sets of all patches
    patches.setSize(nPatches);
    patchTypes.setSize(nPatches);
    patchNames.setSize(nPatches);


    //old polyBoundary
    const polyBoundaryMesh& oPatches = pShapeMesh.boundaryMesh();
    // new patches.
    DynamicList<polyPatch*> newPatches(nPatches);

    // For every boundary face the old patch.
    labelList facePatchID(pShapeMesh.nFaces()-pShapeMesh.nInternalFaces(), -1);
    label cMeshFace = pShapeMesh.nInternalFaces();
    label nBoundaries = 0;


    forAll(patches, patchi)
    {
        const faceList& bFaces = patches[patchi];

        label sz = bFaces.size();
        labelList meshFaces(sz,-1);

        // Search faces by point matching
        forAll(bFaces, j)
        {
            const face& f = bFaces[j];
            label cMeshFace = findFace(pShapeMesh, f);
            meshFaces[j] = cMeshFace;
        }


        //check if patch is internal
        //also check internal/external-ness of first patch face
        //internal faces cannot become foam boundaries
        //if a face is defined as internal but is actually external
        //it will be put in a default wall boundary
        //internal boundaries are simply ignored

        if
        (
            patchTypes[patchi] != "internal"
         && !pShapeMesh.isInternalFace(meshFaces[0])
        )
        {
            //first face is external and has valid non-internal type

            //check all faces for externalness just to be sure
            //and mark patch number to global list
            forAll(meshFaces, i)
            {
                label facei = meshFaces[i];

                if (pShapeMesh.isInternalFace(facei))
                {
                    FatalErrorInFunction
                        << "Face " << facei << " on new patch "
                        << patchNames[patchi]
                        << " is not an external face of the mesh." << endl
                        << exit(FatalError);
                }

                if (facePatchID[facei - pShapeMesh.nInternalFaces()]!= -1)
                {
                    FatalErrorInFunction
                        << "Face " << facei << " on new patch "
                        << patchNames[patchi]
                        << " has already been marked for repatching to"
                        << " patch "
                        << facePatchID[facei - pShapeMesh.nInternalFaces()]
                        << exit(FatalError);
                }
                facePatchID[facei - pShapeMesh.nInternalFaces()] = nBoundaries;
            }

            //add to boundary patch

            Info<< "Adding new patch " << patchNames[patchi]
                << " of type " << patchTypes[patchi]
                << " as patch " << nBoundaries << endl;

            // Add patch to new patch list
            newPatches.append
            (
                polyPatch::New
                (
                    patchTypes[patchi],
                    patchNames[patchi],
                    sz,
                    cMeshFace,
                    nBoundaries,
                    oPatches
                ).ptr()
            );
            nBoundaries++;
            cMeshFace += sz;
        }
        else
        {
            Info<< "Patch " << patchNames[patchi]
                 << " is internal to the mesh "
                 << " and is not being added to the boundary."
                 << endl;

        }
    }

    // Check for any remaining boundary faces
    // and add them to a default wall patch
    // this routine should generally not be invoked
    {
        DynamicList<label> defaultBoundaryFaces(facePatchID.size());
        forAll(facePatchID, idI)
        {
            if (facePatchID[idI] == -1)
            {
                defaultBoundaryFaces.append(idI);
                facePatchID[idI] = nBoundaries;
            }
        }
        defaultBoundaryFaces.shrink();

        if (defaultBoundaryFaces.size())
        {
            Warning << " fluent mesh has " << defaultBoundaryFaces.size()
                    << " undefined boundary faces." << endl
                    << " Adding undefined faces to new patch `default_wall`"
                    << endl;

            // Add patch to new patch list

            newPatches.append
            (
                polyPatch::New
                (
                    wallPolyPatch::typeName,
                    "default_wall",
                    defaultBoundaryFaces.size(),
                    cMeshFace,
                    nBoundaries,
                    oPatches
                ).ptr()
            );
            nBoundaries++;
            cMeshFace += defaultBoundaryFaces.size();
        }
    }

    newPatches.shrink();

    // Use facePatchIDs map to reorder boundary faces into compact regions

    repatchPolyTopoChanger repatcher(pShapeMesh);

    // Add new list of patches
    repatcher.changePatches(newPatches);

    // Change patch ids
    forAll(facePatchID, idI)
    {
        label facei = idI + pShapeMesh.nInternalFaces();

        repatcher.changePatchID(facei, facePatchID[idI]);
    }
    repatcher.repatch();


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));


    // Re-do face matching to write sets
    if (writeSets)
    {
        forAll(patches, patchi)
        {
            const faceList& bFaces = patches[patchi];
            label sz = bFaces.size();

            faceSet pFaceSet(pShapeMesh, patchNames[patchi], sz);

            forAll(bFaces, j)
            {
                const face& f = bFaces[j];
                label cMeshFace = findFace(pShapeMesh, f);
                pFaceSet.insert(cMeshFace);
            }
            Info<< "Writing patch " << patchNames[patchi]
                << " of size " << sz << " to faceSet" << endl;

            pFaceSet.instance() = pShapeMesh.instance();
            pFaceSet.write();
        }
    }


    // Zones
    // will write out cell zones and internal faces for those zones
    // note: zone boundary faces are not added to face zones
    // the names of boundaries bordering on cell zones are written to
    // a list containing the boundary name and cellzone it borders on
    // interior boundaries are handled via faceSets
    // cell zones will only be written if there is more than one

    if (writeZones && cellGroupZoneID.size() > 1)
    {
        Info<< "Adding Zones" << endl;
        List<pointZone*> pz(0);

        label nrCellZones = cellGroupZoneID.size();
        List<cellZone*> cz(nrCellZones);

        // Make face zones for cell zones
        List<faceZone*> fz(nrCellZones);

        // List of patch names and the cellZone(s) they border
        // this is just an info file to make MRF easier to setup
        List<DynamicList<word>> boundaryZones
        (
            pShapeMesh.boundaryMesh().size()
        );

        const polyBoundaryMesh& bPatches = pShapeMesh.boundaryMesh();
        forAll(bPatches, pI)
        {
            boundaryZones[pI].append(bPatches[pI].name());
        }

        label cnt=0;
        SLList<label>::iterator cg = cellGroupZoneID.begin();
        SLList<label>::iterator start = cellGroupStartIndex.begin();
        SLList<label>::iterator end = cellGroupEndIndex.begin();

        for (; cg != cellGroupZoneID.end(); ++cg, ++start, ++end)
        {
            const word& name = patchNameIDs[cg()];
            const word& type = patchTypeIDs[cg()];

            Info<< "Writing cell zone: " << name
                << " of type " << type << " starting at " << start() - 1
                << " ending at " << end() - 1 << " to cellSet." << endl;

            labelList cls(end() - start() + 1);

            // Mark zone cells, used for finding faces
            boolList zoneCell(pShapeMesh.nCells(), false);

            // shift cell indizes by 1
            label nr=0;
            for (label celli = (start() - 1); celli < end(); celli++)
            {
                cls[nr]=celli;
                zoneCell[celli] = true;
                nr++;
            }

            cz[cnt] = new cellZone
            (
                name,
                cls,
                cnt,
                pShapeMesh.cellZones()
            );

            DynamicList<label> zoneFaces(pShapeMesh.nFaces());
            forAll(pShapeMesh.faceNeighbour(), facei)
            {
                label nei = pShapeMesh.faceNeighbour()[facei];
                label own = pShapeMesh.faceOwner()[facei];
                if (nei != -1)
                {
                    if (zoneCell[nei] && zoneCell[own])
                    {
                        zoneFaces.append(facei);
                    }
                }
            }
            zoneFaces.shrink();

            fz[cnt] = new faceZone
            (
                name,
                zoneFaces,
                boolList(zoneFaces.size(), false),
                cnt,
                pShapeMesh.faceZones()
            );

            // Add cell zones to patch zone list
            forAll(bPatches, pI)
            {
                const labelList& faceCells = bPatches[pI].faceCells();
                forAll(faceCells, fcI)
                {
                    if (zoneCell[faceCells[fcI] ])
                    {
                        boundaryZones[pI].append(name);
                        break;
                    }
                }
            }

            cnt++;
        }

        pShapeMesh.addZones(pz, fz, cz);

        forAll(bPatches, pI)
        {
            boundaryZones[pI].shrink();
        }

        fileName bczf
        (
            runTime.path()/runTime.constant()
           /"polyMesh"/"boundaryAdjacentCellZones"
        );

        OFstream boundaryCellZonesFile(bczf);
        forAll(boundaryZones, bzI)
        {
            forAll(boundaryZones[bzI], bzII)
            {
                boundaryCellZonesFile << boundaryZones[bzI][bzII] << " ";
            }

            boundaryCellZonesFile << endl;
        }
    }

    Info<< endl << "Writing mesh..." << flush;

    Info<< " to " << pShapeMesh.instance()/pShapeMesh.meshDir()
        << "  " << flush;

    pShapeMesh.setInstance(pShapeMesh.instance());
    pShapeMesh.write();
    Info<< "done." << endl << endl;

    // Write cellSets for Fluent regions
    // allows easy post-processing
    // set and zone functionality will be integrated some time
    // soon negating the need for double output
    if (writeSets)
    {
        if (cellGroupZoneID.size() > 1)
        {
            Info<< "Writing cell sets" << endl;

            SLList<label>::iterator cg = cellGroupZoneID.begin();
            SLList<label>::iterator start = cellGroupStartIndex.begin();
            SLList<label>::iterator end = cellGroupEndIndex.begin();

            // Note: cellGroupXXX are all Fluent indices (starting at 1)
            // so offset before using.

            for (; cg != cellGroupZoneID.end(); ++cg, ++start, ++end)
            {
                const word& name=patchNameIDs[cg()];
                const word& type=patchTypeIDs[cg()];

                Info<< "Writing cell set: " << name
                    << " of type " << type << " starting at " << start() - 1
                    << " ending at " << end() - 1 << " to cellSet." << endl;

                cellSet internal(pShapeMesh, name, end() - start());

                // shift cell indizes by 1
                for (label celli = start() - 1; celli <= end() - 1; celli++)
                {
                    internal.insert(celli);
                }

                internal.write();
            }
        }
        else
        {
            Info<< "Only one cell group: no set written\n";
        }
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of fluentMeshToFoam.L
 \* ------------------------------------------------------------------------- */

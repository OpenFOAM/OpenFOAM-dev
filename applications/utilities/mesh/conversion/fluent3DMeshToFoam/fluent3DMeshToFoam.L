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
    fluent3DMeshToFoam

Description
    Converts a Fluent mesh to OpenFOAM format.

\*---------------------------------------------------------------------------*/

%{
#undef yyFlexLexer

 /* ------------------------------------------------------------------------ *\
   ------ local definitions
 \* ------------------------------------------------------------------------ */

#include "cyclicPolyPatch.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyMeshZipUpCells.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "oldCyclicPolyPatch.H"
#include "Swap.H"
#include "IFstream.H"
#include "readHexLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Line number
label lineNo = 1;

// Scale factor used to scale points (optional command line argument)
scalar scaleFactor = 1.0;

label dimensionOfGrid = 0;
label nPoints = 0;
label nFaces = 0;
label nCells = 0;

bool hangingNodes = false;

pointField points(0);
faceList faces(0);
labelList owner(0);
labelList neighbour(0);

// Group type and name
Map<word> groupType(100);
Map<word> groupName(100);

// Point groups
DynamicList<label> pointGroupZoneID;
DynamicList<label> pointGroupStartIndex;
DynamicList<label> pointGroupEndIndex;

// Face groups
DynamicList<label> faceGroupZoneID;
DynamicList<label> faceGroupStartIndex;
DynamicList<label> faceGroupEndIndex;

// Cell groups
DynamicList<label> cellGroupZoneID;
DynamicList<label> cellGroupStartIndex;
DynamicList<label> cellGroupEndIndex;
DynamicList<label> cellGroupType;

// Special parsing of (incorrect) Cubit files
bool cubitFile = false;


void uniquify(word& name, HashSet<word>& patchNames)
{
    if (!patchNames.found(name))
    {
        patchNames.insert(name);
    }
    else
    {
        Info<< "    name " << name << " already used";

        label i = 1;
        word baseName = name;

        do
        {
            name = baseName + "_" + Foam::name(i++);
        } while (patchNames.found(name));

        Info<< ", changing to " << name << endl;
    }
}


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
commaPipe                  [,\|]

schemeSpecialInitial       [!$%&*/\\:<=>?~_^#.@']
schemeSpecialSubsequent    [.+-]
schemeSymbol               (({some_space}|{alpha}|{quote}|{schemeSpecialInitial})({alpha}|{quote}|{digit}|{schemeSpecialInitial}|{schemeSpecialSubsequent})*)


identifier                 {alpha}({alpha}|{digit})*
integer                    {decDigit}+
label                      [1-9]{decDigit}*
hexLabel                   {hexDigit}+
zeroLabel                  {digit}*
signedInteger              [-+]?{integer}
word                       ({alpha}|{digit}|{dotColonDash})*
wordBraces                 ({word}|{lbrac}|{rbrac})*
wordBracesExtras           ({word}|{lbrac}|{rbrac}|{commaPipe})*

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
textBraces                 ({space}({wordBraces}*{space})*)
textExtras                 ({space}({word}*{commaPipe}{space})*)
textBracesExtras           ({space}({wordBracesExtras}*{space})*)
anythingInBlock            ([^)]*)

dateDDMMYYYY               ({digit}{digit}"/"{digit}{digit}"/"{digit}{digit}{digit}{digit})
dateDDMonYYYY              ((({digit}{digit}{space})|({digit}{space})){alpha}*{space}{digit}{digit}{digit}{digit})
time                       ({digit}{digit}":"{digit}{digit}":"{digit}{digit})

versionNumber              ({digit}|".")*

header                     {space}"(1"{space}
dimension                  {space}"(2"{space}
points                     {space}"(10"{space}
faces                      {space}"(13"{space}
cells                      {space}"(12"{space}
zoneVariant1               {space}"(39"{space}
zoneVariant2               {space}"(45"{space}
faceTree                   {space}"(59"{space}

comment                    "0"{space}
unknownPeriodicFace        "17"{space}
periodicFace               "18"{space}
cellTree                   "58"{space}
faceParents                "61"{space}
ignoreBlocks               ("4"|"37"|"38"|"40"|"41"|"60"|"64"){space}

redundantBlock             {space}({comment}|{unknownPeriodicFace}|{periodicFace}|{cellTree}|{faceParents}|{ignoreBlocks}){space}

endOfSection               {space}")"{space}



 /* ------------------------------------------------------------------------ *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------ */

%option stack

%x readHeader
%x readDimension
%x readPoint
%x readPointHeader
%x readNumberOfPoints
%x readPointGroupData
%x readPointData
%x readScalarList
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
%x readCellsUniform
%x zone
%x readZoneHeader
%x readZoneGroupData
%x readZoneData
%x readZoneBlock

%x ignoreBlock
%x ignoreEmbeddedBlock
%%

%{
    // End of read character pointer returned by strtol and strtod
    char* endPtr;

    // Point data
    label pointGroupNumberOfComponents = 3;
    label pointi = 0; // index used for reading points
    label cmpt = 0;   // component index used for reading points

    // Face data
    label faceGroupElementType = -1;
    label facei = 0;  // index used for reading faces
%}


 /* ------------------------------------------------------------------------ *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------ */

 /*                      ------ Reading control header ------                */

{header} {
        BEGIN(readHeader);
    }

<readHeader>{quote}{textBracesExtras}{quote} {
        Info<< "Reading header: " << YYText() << endl;
    }

{dimension} {
        BEGIN(readDimension);
    }

<readDimension>{space}{label}{space} {
        dimensionOfGrid = atoi(YYText());
        Info<< "Dimension of grid: " << dimensionOfGrid << endl;
    }


{points} {
        BEGIN(readPointHeader);
    }

<readPointHeader>{space}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfPoints);
    }

<readNumberOfPoints>{hexLabel}{space}{labelList} {
        nPoints = strtol(YYText(), &endPtr, 16);
        Info<< "Number of points: " << nPoints << endl;
        points.setSize(nPoints);

        // Ignore rest of stream
    }

<readPointHeader>{space}{lbrac} {
        BEGIN(readPointGroupData);
    }

<readPointGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{labelList} {
        // Read point zone-ID, start and end-label
        // the indices will be used for checking later.
        pointGroupZoneID.append(strtol(YYText(), &endPtr, 16));

        // In FOAM, indices start from zero - adjust
        pointGroupStartIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        pointGroupEndIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        // point group type skipped
        (void)strtol(endPtr, &endPtr, 16);

        pointi = pointGroupStartIndex.last();

        // reset number of components to default
        pointGroupNumberOfComponents = 3;

        // read number of components in the vector
        if (endPtr < &(YYText()[YYLeng()-1]))
        {
            pointGroupNumberOfComponents = strtol(endPtr, &endPtr, 16);
        }

        Info<< "PointGroup: "
            << pointGroupZoneID.last()
            << " start: "
            << pointGroupStartIndex.last()
            << " end: "
            << pointGroupEndIndex.last() << flush;
    }

<readNumberOfPoints,readPointGroupData>{endOfSection} {
        BEGIN(readPointData);
    }

<readPointData>{space}{lbrac}{space} {
        Info<< ".  Reading points..." << flush;
        cmpt = 0;
        yy_push_state(readScalarList);
    }

<readScalarList>{signedInteger}{space} {
        points[pointi][cmpt++] = scaleFactor*atol(YYText());

        if (cmpt == pointGroupNumberOfComponents)
        {
            if (pointGroupNumberOfComponents == 2)
            {
                points[pointi].z() = 0.0;
            }

            cmpt = 0;
            pointi++;
        }
    }

<readScalarList>{scalar}{space} {
        points[pointi][cmpt++] = scaleFactor*atof(YYText());

        if (cmpt == pointGroupNumberOfComponents)
        {
            if (pointGroupNumberOfComponents == 2)
            {
                points[pointi].z() = 0.0;
            }

            cmpt = 0;
            pointi++;
        }
    }

<readScalarList>{endOfSection} {
        Info<< "done." << endl;

        // check read of points
        if (pointi != pointGroupEndIndex.last()+1)
        {
            Warning
                << "Problem with reading points: " << nl
                << "    start index: "
                << pointGroupStartIndex.last()
                << " end index: "
                << pointGroupEndIndex.last()
                << " last points read: " << pointi << nl
                << "    on line " << lineNo << endl;
        }

        yy_pop_state();
    }

{faces} {
        BEGIN(readFaceHeader);
    }

<readFaceHeader>{space}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfFaces);
    }

<readNumberOfFaces>{space}{hexLabel}{space}{labelListElement}+ {
        nFaces = strtol(YYText(), &endPtr, 16);

        Info<< "Number of faces: " << nFaces << endl;

        faces.setSize(nFaces);
        owner.setSize(nFaces);
        neighbour.setSize(nFaces);

        // Type and element type not read
    }

<readFaceHeader>{space}{lbrac} {
        BEGIN(readFaceGroupData);
    }

<readFaceGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{hexLabelListElement}+ {
        // read fluentFace zone-ID, start and end-label
        faceGroupZoneID.append(strtol(YYText(), &endPtr, 16));

        // In FOAM, indices start from zero - adjust
        faceGroupStartIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        faceGroupEndIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        // face group type
        (void)strtol(endPtr, &endPtr, 16);

        faceGroupElementType = strtol(endPtr, &endPtr, 16);

        facei = faceGroupStartIndex.last();

        Info<< "FaceGroup: "
            << faceGroupZoneID.last()
            << " start: "
            << faceGroupStartIndex.last()
            << " end: "
            << faceGroupEndIndex.last() << flush;
    }

<readNumberOfFaces,readFaceGroupData>{space}{endOfSection} {
        BEGIN(readFaceData);
    }

<readFaceData>{space}{lbrac} {
        if (faceGroupElementType == 0 || faceGroupElementType > 4)
        {
            Info<< ".  Reading mixed faces..." << flush;
            yy_push_state(readFacesMixed);
        }
        else
        {
            Info<< ".  Reading uniform faces..." << flush;
            yy_push_state(readFacesUniform);
        }
    }

<readFacesMixed>{space}{hexLabelList} {
        face& curFaceLabels = faces[facei];

        // set size of label list
        curFaceLabels.setSize(strtol(YYText(), &endPtr, 16));

        forAll(curFaceLabels, i)
        {
            curFaceLabels[i] = strtol(endPtr, &endPtr, 16) - 1;
        }

        // read neighbour and owner. Neighbour comes first
        neighbour[facei] = strtol(endPtr, &endPtr, 16) - 1;
        owner[facei] = strtol(endPtr, &endPtr, 16) - 1;
        facei++;
    }

<readFacesUniform>{space}{hexLabelList} {
        face& curFaceLabels = faces[facei];

        // Set size of label list.
        curFaceLabels.setSize(faceGroupElementType);

        curFaceLabels[0] = strtol(YYText(), &endPtr, 16) - 1;

        for (int i=1; i<faceGroupElementType; i++)
        {
            curFaceLabels[i] = strtol(endPtr, &endPtr, 16) - 1;
        }

        // read neighbour and owner. Neighbour comes first
        neighbour[facei] = strtol(endPtr, &endPtr, 16) - 1;
        owner[facei] = strtol(endPtr, &endPtr, 16) - 1;
        facei++;
    }

<readFacesMixed,readFacesUniform>{space}{endOfSection} {
        Info<< "done." << endl;

        // check read of fluentFaces
        if (facei != faceGroupEndIndex.last()+1)
        {
            Warning
                << "Problem with reading fluentFaces: " << nl
                << "    start index: "
                << faceGroupStartIndex.last()
                << " end index: "
                << faceGroupEndIndex.last()
                << " last fluentFaces read: " << facei << nl
                << "    on line " << lineNo << endl;
        }

        yy_pop_state();
    }


{cells} {
        BEGIN(readCellHeader);
    }

<readCellHeader>{space}{lbrac}{space}"0"{space}"1"{space} {
        BEGIN(readNumberOfCells);
    }

<readNumberOfCells>{space}{hexLabel}{space}{labelListElement}+ {
        nCells = strtol(YYText(), &endPtr, 16);
        Info<< "Number of cells: " << nCells << endl;
    }

<readCellHeader>{space}{lbrac} {
        BEGIN(readCellGroupData);
    }

<readCellGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel} {
        // Warning. This entry must be above the next one because of the lexing
        // rules. It is introduced to deal with the problem of reading
        // non-standard cell definition from Tgrid, which misses the type label.

        Warning
            << "Tgrid syntax problem: " << YYText() << nl
            << "    on line " << lineNo << endl;

        // read cell zone-ID, start and end-label
        cellGroupZoneID.append(strtol(YYText(), &endPtr, 16));

        // the indices will be used for checking later.
        cellGroupStartIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        cellGroupEndIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        cellGroupType.append(strtol(endPtr, &endPtr, 16));

        Info<< "CellGroup: "
            << cellGroupZoneID.last()
            << " start: "
            << cellGroupStartIndex.last()
            << " end: "
            << cellGroupEndIndex.last()
            << " type: "
            << cellGroupType.last()
            << endl;
    }

<readCellGroupData>{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel}{space}{hexLabel} {
        // Warning. See above

        // read cell zone-ID, start and end-label
        cellGroupZoneID.append(strtol(YYText(), &endPtr, 16));

        // the indices will be used for checking later.
        cellGroupStartIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        cellGroupEndIndex.append(strtol(endPtr, &endPtr, 16) - 1);

        cellGroupType.append(strtol(endPtr, &endPtr, 16));

        // Note. Potentially skip cell set if type is zero.
        (void)strtol(endPtr, &endPtr, 16);

        Info<< "CellGroup: "
            << cellGroupZoneID.last()
            << " start: "
            << cellGroupStartIndex.last()
            << " end: "
            << cellGroupEndIndex.last()
            << " type: "
            << cellGroupType.last()
            << endl;
    }

<readNumberOfCells,readCellGroupData>{endOfSection} {
        BEGIN(readCellData);
    }

<readCellData>{space}{lbrac} {
        // Quickly scan to the end of the cell data block and discard
        int c;
        while ((c = yyinput()) != 0 && c != ')')
        {}
    }

{faceTree} {
        // There are hanging nodes in the mesh so make sure it gets zipped-up
        hangingNodes = true;
        yy_push_state(ignoreBlock);
    }

{zoneVariant1} {
        BEGIN(readZoneHeader);
    }

{zoneVariant2} {
        BEGIN(readZoneHeader);
    }

<readZoneHeader>{space}{lbrac} {
        BEGIN(readZoneGroupData);
    }

<readZoneGroupData>{space}{hexLabel}{space}{word}{space}{word}{space}{label}? {
        IStringStream zoneDataStream(YYText());

        // cell zone-ID not in hexadecimal!!! Inconsistency
        label zoneID = -1;

        if (cubitFile)
        {
            zoneID = readHexLabel(zoneDataStream);
        }
        else
        {
            zoneID = readLabel(zoneDataStream);
        }

        groupType.insert(zoneID, word(zoneDataStream));
        groupName.insert(zoneID, word(zoneDataStream));

        Info<< "Zone: " << zoneID
            << " name: " << groupName[zoneID]
            << " type: " << groupType[zoneID] << flush;
    }

<readZoneGroupData>{endOfSection} {
        BEGIN(readZoneData);
    }

<readZoneData>{space}{lbrac} {
        Info<< ".  Reading zone data..." << flush;
        yy_push_state(readZoneBlock);
    }

<readZoneBlock>{space}{schemeSymbolList} {
    }

<readZoneBlock>{lbrac} {
        // Warning
        //    << "Found unknown block in zone: " << YYText() << nl
        //    << "    on line " << lineNo << endl;
        yy_push_state(ignoreBlock);
    }

<readZoneBlock>{endOfSection} {
        Info<< "done." << endl;
        yy_pop_state();
    }



 /*             ------ Reading end of section and others ------              */

<readHeader,readDimension,readPointData,readFaceData,readCellData,readZoneData>{space}{endOfSection} {
        BEGIN(INITIAL);
    }

 /*    ------ Reading unknown type or non-standard comment ------            */

{lbrac}{label} {
        Warning
            << "Found unknown block of type: "
            << Foam::string(YYText())(1, YYLeng()-1) << nl
            << "    on line " << lineNo << endl;

        yy_push_state(ignoreBlock);
    }

{lbrac}{redundantBlock} {
        yy_push_state(ignoreBlock);
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{quote}{text}{quote} {
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{schemeSymbol} {
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{lbrac} {
        yy_push_state(ignoreEmbeddedBlock);

    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{endOfSection} {
        yy_pop_state();
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{labelList} {
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{hexLabelList} {
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{scalarList} {
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{schemeSymbolList} {
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{text} {
    }

<ignoreBlock,ignoreEmbeddedBlock>{space}{textExtras} {
    }

 /* ------              Count newlines.                              ------  */

<*>\n {
        lineNo++;
    }


 /* ------              Ignore remaining space.                      ------  */

<*>{some_space}|\r {
    }


 /* ------              Any other characters are errors.              ------ */

<*>. {
        // This is a catch all.
        FatalErrorInFunction
            << "Do not understand characters: " << YYText() << nl
            << "    on line " << lineNo
            << exit(FatalError);
    }


 /*  ------ On EOF return to previous file, if none exists terminate. ------ */

<<EOF>> {
            yyterminate();
    }
%%

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
        "ignoreCellGroups",
        "names",
        "specify cell groups to ignore"
    );
    argList::addOption
    (
        "ignoreFaceGroups",
        "names",
        "specify face groups to ignore"
    );

    argList::addBoolOption
    (
        "cubit",
        "special parsing of (incorrect) cubit files"
    );

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    args.optionReadIfPresent("scale", scaleFactor);

    wordHashSet ignoreCellGroups;
    wordHashSet ignoreFaceGroups;

    args.optionReadIfPresent("ignoreCellGroups", ignoreCellGroups);
    args.optionReadIfPresent("ignoreFaceGroups", ignoreFaceGroups);

    cubitFile = args.options().found("cubit");

    if (cubitFile)
    {
        Info<< nl
            << "Assuming Cubit generated file"
            << " (incorrect face orientation; hexadecimal zoneIDs)."
            << nl << endl;
    }


    #include "createTime.H"

    const fileName fluentFile = args[1];
    IFstream fluentStream(fluentFile);

    if (!fluentStream)
    {
        FatalErrorInFunction
            << ": file " << fluentFile << " not found"
            << exit(FatalError);
    }

    yyFlexLexer lexer(&fluentStream.stdStream());

    while (lexer.yylex() != 0)
    {}

    Info<< "\nFINISHED LEXING\n\n";

    if (dimensionOfGrid != 3)
    {
        FatalErrorInFunction
            << "Mesh is not 3D, dimension of grid: " << dimensionOfGrid
            << exit(FatalError);
    }

    pointGroupZoneID.shrink();
    pointGroupStartIndex.shrink();
    pointGroupEndIndex.shrink();

    faceGroupZoneID.shrink();
    faceGroupStartIndex.shrink();
    faceGroupEndIndex.shrink();

    cellGroupZoneID.shrink();
    cellGroupStartIndex.shrink();
    cellGroupEndIndex.shrink();
    cellGroupType.shrink();


    // Pre-filtering: flip "owner" boundary or wrong oriented internal
    // faces and move to neighbour

    boolList fm(faces.size(), false);
    forAll(faces, facei)
    {
        if
        (
            owner[facei] == -1
         || (neighbour[facei] != -1 && owner[facei] > neighbour[facei])
        )
        {
            fm[facei] = true;
            if (!cubitFile)
            {
                faces[facei].flip();
            }
            Swap(owner[facei], neighbour[facei]);
        }
    }


    // Foam type for Fluent type
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    HashTable<word> fluentToFoamType;

    fluentToFoamType.insert("pressure", polyPatch::typeName);
    fluentToFoamType.insert("pressure-inlet", polyPatch::typeName);
    fluentToFoamType.insert("inlet-vent", polyPatch::typeName);
    fluentToFoamType.insert("intake-fan", polyPatch::typeName);
    fluentToFoamType.insert("pressure-outlet", polyPatch::typeName);
    fluentToFoamType.insert("exhaust-fan", polyPatch::typeName);
    fluentToFoamType.insert("outlet-vent", polyPatch::typeName);
    fluentToFoamType.insert("pressure-far-field", polyPatch::typeName);
    fluentToFoamType.insert("velocity-inlet", polyPatch::typeName);
    fluentToFoamType.insert("mass-flow-inlet", polyPatch::typeName);
    fluentToFoamType.insert("outflow", polyPatch::typeName);

    fluentToFoamType.insert("wall" , wallPolyPatch::typeName);

    fluentToFoamType.insert("symmetry",  symmetryPolyPatch::typeName);
    fluentToFoamType.insert("axis",  symmetryPolyPatch::typeName);

    fluentToFoamType.insert("interior", polyPatch::typeName);
    fluentToFoamType.insert("interface", polyPatch::typeName);
    fluentToFoamType.insert("internal", polyPatch::typeName);
    fluentToFoamType.insert("solid", polyPatch::typeName);
    fluentToFoamType.insert("fan", oldCyclicPolyPatch::typeName);
    fluentToFoamType.insert("radiator", polyPatch::typeName);
    fluentToFoamType.insert("porous-jump", polyPatch::typeName);

    //- Periodic halves map directly into split cyclics. The problem is the
    //  initial matching since we require knowledge of the transformation.
    //  It is ok if the periodics are already ordered. We should read the
    //  periodic shadow faces section (section 18) to give use the ordering
    //  For now just disable.
    // fluentToFoamType.insert("periodic", cyclicPolyPatch::typeName);


    // Foam patch type for Fluent zone type
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HashSet<word> fluentGroupToFoamPatch;
    fluentGroupToFoamPatch.insert("wall");
    fluentGroupToFoamPatch.insert("fan");


    // Create initial empty polyMesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferCopy(pointField()),
        xferCopy(faceList()),
        xferCopy(labelList()),
        xferCopy(labelList())
    );


    // Check the cell groups for zones ignoring those in ignoreCellGroups
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    label nCellZones = 0;
    labelList cellZoneIDs(cellGroupZoneID.size());

    forAll(cellGroupZoneID, cgi)
    {
        if (!ignoreCellGroups.found(groupName[cellGroupZoneID[cgi] ]))
        {
            cellZoneIDs[nCellZones++] = cgi;
        }
    }

    cellZoneIDs.setSize(nCellZones);


    // Check the face groups for boundary patches, baffles and faceZones
    // ignoring the interior zones in ignoreCellGroups
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DynamicList<label> patchIDs(faceGroupZoneID.size());
    DynamicList<label> faceZoneIDs(faceGroupZoneID.size());

    forAll(faceGroupZoneID, fgi)
    {
        label zoneID = faceGroupZoneID[fgi];
        label start = faceGroupStartIndex[fgi];

        if (groupType.found(zoneID))
        {
            const word& type = groupType[zoneID];

            // Check the first element of neighbour for boundary group
            if (neighbour[start] == -1 || fluentGroupToFoamPatch.found(type))
            {
                patchIDs.append(fgi);
            }
            else
            {
                if (!ignoreFaceGroups.found(groupName[faceGroupZoneID[fgi] ]))
                {
                    faceZoneIDs.append(fgi);
                }
            }
        }
        else if (hangingNodes)
        {
            label end = faceGroupEndIndex[fgi];

            Info<< "Unknown FaceGroup " << zoneID
                << " assumed to be parent faces of refinement "
                   "patterns and ignored."
                << endl;

            // Set the owner of these faces to -1 so that they do not get
            // added to the mesh
            for (label facei = start; facei <= end; facei++)
            {
                owner[facei] = -1;
            }
        }
        else
        {
            if (neighbour[start] == -1)
            {
                // Boundary face in unknown group. Create a patch for it.
                groupType.insert(zoneID, "unknown");
                groupName.insert(zoneID, "FaceGroup" + Foam::name(zoneID));
                patchIDs.append(fgi);
                Info<< "Created patch " << fgi << " for unknown FaceGroup "
                    << zoneID << '.' << endl;
            }
            else
            {
                WarningInFunction
                    << "Unknown FaceGroup " << zoneID << " not in a zone"
                    << endl;
            }
        }
    }

    patchIDs.shrink();
    faceZoneIDs.shrink();


    // Add empty patches
    // ~~~~~~~~~~~~~~~~~

    List<polyPatch*> newPatches(patchIDs.size());
    HashSet<word> patchNames;

    forAll(patchIDs, patchi)
    {
        label zoneID = faceGroupZoneID[patchIDs[patchi] ];
        word name = groupName[zoneID];
        const word& type = groupType[zoneID];

        Info<< "Creating patch " << patchi
            << " for zone: " << zoneID
            << " name: " << name
            << " type: " << type
            << endl;

        uniquify(name, patchNames);

        HashTable<word>::const_iterator iter = fluentToFoamType.find(type);

        if (iter != fluentToFoamType.end())
        {
            // See if we have a periodic and can derive the other side.
            word neighbPatchName;
            if (iter() == cyclicPolyPatch::typeName)
            {
                // Periodic
                size_t n = name.rfind("-SIDE-1");

                if (n != string::npos)
                {
                    neighbPatchName = name.substr(0, n) + "-SIDE-2";
                }
                else
                {
                    n = name.rfind("-SIDE-2");
                    if (n != string::npos)
                    {
                        neighbPatchName = name.substr(0, n) + "-SIDE-1";
                    }
                }
            }

            if (neighbPatchName.size())
            {
                Info<< "Adding cyclicPolyPatch for Fluent zone " << name
                    << " with neighbour patch " << neighbPatchName
                    << endl;

                newPatches[patchi] = new cyclicPolyPatch
                (
                    name,
                    0,
                    0,
                    patchi,
                    mesh.boundaryMesh(),
                    neighbPatchName,
                    cyclicPolyPatch::NOORDERING,
                    Zero,
                    Zero,
                    Zero
                );
            }
            else
            {
                newPatches[patchi] = polyPatch::New
                (
                    iter(),
                    name,
                    0,
                    0,
                    patchi,
                    mesh.boundaryMesh()
                ).ptr();
            }
        }
        else
        {
            Info<< "Adding polyPatch for unknown Fluent type " << type
                << endl;

            newPatches[patchi] = new polyPatch
            (
                name,
                0,
                0,
                patchi,
                mesh.boundaryMesh(),
                polyPatch::typeName
            );
        }
    }
    mesh.addPatches(newPatches);


    // Add empty zones
    // ~~~~~~~~~~~~~~~

    // Cell zones
    mesh.cellZones().setSize(cellZoneIDs.size());
    HashSet<word> cellZoneNames;

    forAll(cellZoneIDs, cellZonei)
    {
        label zoneID = cellGroupZoneID[cellZoneIDs[cellZonei] ];
        word name = groupName[zoneID];
        const word& type = groupType[zoneID];

        Info<< "Creating cellZone " << cellZonei
            << " name: " << name
            << " type: " << type
            << endl;

        uniquify(name, cellZoneNames);

        mesh.cellZones().set
        (
            cellZonei,
            new cellZone
            (
                name,
                labelList(0),
                cellZonei,
                mesh.cellZones()
            )
        );
    }

    // Face zones
    mesh.faceZones().setSize(faceZoneIDs.size());
    HashSet<word> faceZoneNames;

    forAll(faceZoneIDs, faceZonei)
    {
        label zoneID = faceGroupZoneID[faceZoneIDs[faceZonei] ];
        word name = groupName[zoneID];
        const word& type = groupType[zoneID];

        Info<< "Creating faceZone " << faceZonei
            << " name: " << name
            << " type: " << type
            << endl;

        uniquify(name, faceZoneNames);

        mesh.faceZones().set
        (
            faceZonei,
            new faceZone
            (
                name,
                labelList(0),
                boolList(0),
                faceZonei,
                mesh.faceZones()
            )
        );
    }


    // Modify mesh for points/cells/faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mesh-change container
    polyTopoChange meshMod(mesh, false);

    // Add all points
    forAll(points, pointi)
    {
        meshMod.addPoint(points[pointi], pointi, -1, true);
    }
    points.setSize(0);

    // Add all cells
    for (label celli = 0; celli < nCells; celli++)
    {
        meshMod.addCell
        (
            -1,         // masterPointID
            -1,         // masterEdgeID
            -1,         // masterFaceID
            celli,      // masterCellID
            -1          // zoneID
        );
    }

    // Modify cells to be in zones as required
    forAll(cellZoneIDs, cellZonei)
    {
        label cgi = cellZoneIDs[cellZonei];

        for
        (
            label celli = cellGroupStartIndex[cgi];
            celli <= cellGroupEndIndex[cgi];
            celli++
        )
        {
            meshMod.modifyCell(celli, cellZonei);
        }
    }


    bool doneWarning = false;

    // Add faceZone faces
    forAll(faceZoneIDs, faceZonei)
    {
        label fgi = faceZoneIDs[faceZonei];
        label start = faceGroupStartIndex[fgi];
        label end = faceGroupEndIndex[fgi];
        label zoneID = faceGroupZoneID[fgi];

        Info<< "faceZone from Fluent indices: " << start
            << " to: " << end
            << " type: " << groupType[zoneID]
            << endl;

        for (label facei = start; facei <= end; facei++)
        {
            if (owner[facei] >= nCells || neighbour[facei] >= nCells)
            {
                if (!doneWarning)
                {
                    WarningInFunction
                        << "Ignoring internal face " << facei
                        << " on FaceZone " << zoneID
                        << " since owner " << owner[facei] << " or neighbour "
                        << neighbour[facei] << " outside range of cells 0.."
                        << nCells-1 << endl
                        << "    Suppressing future warnings." << endl;
                    doneWarning = true;
                }
            }
            else
            {
                meshMod.addFace
                (
                    faces[facei],
                    owner[facei],
                    neighbour[facei],
                    -1,                 // masterPointID
                    -1,                 // masterEdgeID
                    facei,              // masterFace
                    false,              // flipFaceFlux
                    -1,                 // patchID
                    faceZonei,          // zoneID
                    fm[facei]           // zoneFlip
                );
            }

            // Mark face as being done
            owner[facei] = -1;
        }
    }

    // Add patch faces
    forAll(patchIDs, patchi)
    {
        label fgi = patchIDs[patchi];
        label start = faceGroupStartIndex[fgi];
        label end = faceGroupEndIndex[fgi];
        label zoneID = faceGroupZoneID[fgi];

        Info<< "patch " << patchi << " from Fluent indices: " << start
            << " to: " << end
            << " type: " << groupType[zoneID]
            << endl;

        for (label facei = start; facei <= end; facei++)
        {
            if (owner[facei] >= nCells || neighbour[facei] >= nCells)
            {
                if (!doneWarning)
                {
                    WarningInFunction
                        << "Ignoring patch face " << facei
                        << " on FaceZone " << zoneID
                        << " since owner " << owner[facei] << " or neighbour "
                        << neighbour[facei] << " outside range of cells 0.."
                        << nCells-1 << endl
                        << "    Suppressing future warnings." << endl;
                    doneWarning = true;
                }
            }
            else
            {
                meshMod.addFace
                (
                    faces[facei],
                    owner[facei],
                    -1,
                    -1,                 // masterPointID
                    -1,                 // masterEdgeID
                    facei,              // masterFace
                    false,              // flipFaceFlux
                    patchi,             // patchID
                    -1,                 // zoneID
                    false               // zoneFlip
                );

                // For baffles create the opposite face
                if (neighbour[start] != -1)
                {
                    meshMod.addFace
                    (
                        faces[facei].reverseFace(),
                        neighbour[facei],
                        -1,
                        -1,                 // masterPointID
                        -1,                 // masterEdgeID
                        facei,              // masterFace
                        false,              // flipFaceFlux
                        patchi,             // patchID
                        -1,                 // zoneID
                        false               // zoneFlip
                    );
                }
            }
            // Mark face as being done
            owner[facei] = -1;
        }
    }

    // Add remaining internal faces
    forAll(owner, facei)
    {
        if (owner[facei] != -1)
        {
            // Check the face being added as an internal face actually is one
            if (neighbour[facei] == -1)
            {
                FatalErrorInFunction
                    << "Attempt of add internal face " << facei
                    << " which is a boundary face"
                    << exit(FatalError);
            }

            if (owner[facei] >= nCells || neighbour[facei] >= nCells)
            {
                if (!doneWarning)
                {
                    WarningInFunction
                        << "Ignoring internal face " << facei
                        << " since owner " << owner[facei] << " or neighbour "
                        << neighbour[facei] << " outside range of cells 0.."
                        << nCells-1 << endl
                        << "    Suppressing future warnings." << endl;
                    doneWarning = true;
                }
            }
            else
            {
                meshMod.addFace
                (
                    faces[facei],
                    owner[facei],
                    neighbour[facei],
                    -1,                 // masterPointID
                    -1,                 // masterEdgeID
                    facei,              // masterFace
                    false,              // flipFaceFlux
                    -1,                 // patchID
                    -1,                 // zoneID
                    false               // zoneFlip
                );
            }
        }
    }

    // Reclaim storage
    faces.setSize(0);
    owner.setSize(0);
    neighbour.setSize(0);


    // Modify mesh
    // ~~~~~~~~~~~

    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Zip-up the mesh if it contained hanging nodes
    if (hangingNodes)
    {
        Info<< "Zipping mesh to remove hanging nodes" << endl;
        polyMeshZipUpCells(mesh);
    }

    mesh.setInstance(runTime.constant());

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing mesh to " << mesh.objectPath() << endl;
    mesh.write();


    Info<< "\nEnd\n" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------ *\
    ------ End of fluentMeshToFoam.L
 \* ------------------------------------------------------------------------ */

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "vtkUnstructuredReader.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "stringIOList.H"
#include "cellModeller.H"
#include "vectorIOField.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(vtkUnstructuredReader, 1);   //0);

    template<>
    const char*
    NamedEnum<vtkUnstructuredReader::vtkDataType, 8>::names[] =
    {
        "int",
        "unsigned_int",
        "long",
        "unsigned_long",
        "float",
        "double",
        "string",
        "vtkIdType"
    };
    const NamedEnum<vtkUnstructuredReader::vtkDataType, 8>
    vtkUnstructuredReader::vtkDataTypeNames;


    template<>
    const char*
    NamedEnum<vtkUnstructuredReader::vtkDataSetType, 3>::names[] =
    {
        "FIELD",
        "SCALARS",
        "VECTORS"
    };
    const NamedEnum<vtkUnstructuredReader::vtkDataSetType, 3>
    vtkUnstructuredReader::vtkDataSetTypeNames;


    template<>
    const char*
    NamedEnum<vtkUnstructuredReader::parseMode, 5>::names[] =
    {
        "NOMODE",
        "UNSTRUCTURED_GRID",
        "POLYDATA",
        "CELL_DATA",
        "POINT_DATA"
    };
    const NamedEnum<vtkUnstructuredReader::parseMode, 5>
    vtkUnstructuredReader::parseModeNames;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkUnstructuredReader::warnUnhandledType
(
    Istream& inFile,
    const label type,
    labelHashSet& warningGiven
) const
{
    if (warningGiven.insert(type))
    {
        IOWarningIn("vtkUnstructuredReader::warnUnhandledType(..)", inFile)
            << "Skipping unknown cell type " << type << endl;
    }
}


// Split cellTypes into cells, faces and lines
void Foam::vtkUnstructuredReader::extractCells
(
    Istream& inFile,
    const labelList& cellTypes,
    const labelList& cellVertData
)
{
    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    labelList tetPoints(4);
    labelList pyrPoints(5);
    labelList prismPoints(6);
    labelList hexPoints(8);

    label cellI = cells_.size();
    cells_.setSize(cellI+cellTypes.size());
    cellMap_.setSize(cells_.size(), -1);

    label faceI = faces_.size();
    faces_.setSize(faceI+cellTypes.size());
    faceMap_.setSize(faces_.size(), -1);

    label lineI = lines_.size();
    lines_.setSize(lineI+cellTypes.size());
    lineMap_.setSize(lines_.size(), -1);

    label dataIndex = 0;


    // To mark whether unhandled type has been visited.
    labelHashSet warningGiven;

    forAll(cellTypes, i)
    {
        switch (cellTypes[i])
        {
            case VTK_VERTEX:
            {
                warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 1)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 1 for VTK_VERTEX but found "
                        << nRead << exit(FatalIOError);
                }
                dataIndex += nRead;
            }
            break;

            case VTK_POLY_VERTEX:
            {
                warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                dataIndex += nRead;
            }
            break;

            case VTK_LINE:
            {
                //warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 2)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 2 for VTK_LINE but found "
                        << nRead << exit(FatalIOError);
                }
                lineMap_[lineI] = i;
                labelList& segment = lines_[lineI++];
                segment.setSize(2);
                segment[0] = cellVertData[dataIndex++];
                segment[1] = cellVertData[dataIndex++];
            }
            break;

            case VTK_POLY_LINE:
            {
                //warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                lineMap_[lineI] = i;
                labelList& segment = lines_[lineI++];
                segment.setSize(nRead);
                forAll(segment, i)
                {
                    segment[i] = cellVertData[dataIndex++];
                }
            }
            break;

            case VTK_TRIANGLE:
            {
                faceMap_[faceI] = i;
                face& f = faces_[faceI++];
                f.setSize(3);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 3)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 3 for VTK_TRIANGLE but found "
                        << nRead << exit(FatalIOError);
                }
                f[0] = cellVertData[dataIndex++];
                f[1] = cellVertData[dataIndex++];
                f[2] = cellVertData[dataIndex++];
            }
            break;

            case VTK_QUAD:
            {
                faceMap_[faceI] = i;
                face& f = faces_[faceI++];
                f.setSize(4);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 4)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 4 for VTK_QUAD but found "
                        << nRead << exit(FatalIOError);
                }
                f[0] = cellVertData[dataIndex++];
                f[1] = cellVertData[dataIndex++];
                f[2] = cellVertData[dataIndex++];
                f[3] = cellVertData[dataIndex++];
            }
            break;

            case VTK_POLYGON:
            {
                faceMap_[faceI] = i;
                face& f = faces_[faceI++];
                label nRead = cellVertData[dataIndex++];
                f.setSize(nRead);
                forAll(f, fp)
                {
                    f[fp] = cellVertData[dataIndex++];
                }
            }
            break;

            case VTK_TETRA:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 4)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 4 for VTK_TETRA but found "
                        << nRead << exit(FatalIOError);
                }
                tetPoints[0] = cellVertData[dataIndex++];
                tetPoints[1] = cellVertData[dataIndex++];
                tetPoints[2] = cellVertData[dataIndex++];
                tetPoints[3] = cellVertData[dataIndex++];
                cellMap_[cellI] = i;
                cells_[cellI++] = cellShape(tet, tetPoints, true);
            }
            break;

            case VTK_PYRAMID:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 5)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 5 for VTK_PYRAMID but found "
                        << nRead << exit(FatalIOError);
                }
                pyrPoints[0] = cellVertData[dataIndex++];
                pyrPoints[1] = cellVertData[dataIndex++];
                pyrPoints[2] = cellVertData[dataIndex++];
                pyrPoints[3] = cellVertData[dataIndex++];
                pyrPoints[4] = cellVertData[dataIndex++];
                cellMap_[cellI] = i;
                cells_[cellI++] = cellShape(pyr, pyrPoints, true);
            }
            break;

            case VTK_WEDGE:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 6)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 6 for VTK_WEDGE but found "
                        << nRead << exit(FatalIOError);
                }
                prismPoints[0] = cellVertData[dataIndex++];
                prismPoints[1] = cellVertData[dataIndex++];
                prismPoints[2] = cellVertData[dataIndex++];
                prismPoints[3] = cellVertData[dataIndex++];
                prismPoints[4] = cellVertData[dataIndex++];
                prismPoints[5] = cellVertData[dataIndex++];
                cellMap_[cellI] = i;
                cells_[cellI++] = cellShape(prism, prismPoints, true);
            }
            break;

            case VTK_HEXAHEDRON:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 8)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 8 for VTK_HEXAHEDRON but found "
                        << nRead << exit(FatalIOError);
                }
                hexPoints[0] = cellVertData[dataIndex++];
                hexPoints[1] = cellVertData[dataIndex++];
                hexPoints[2] = cellVertData[dataIndex++];
                hexPoints[3] = cellVertData[dataIndex++];
                hexPoints[4] = cellVertData[dataIndex++];
                hexPoints[5] = cellVertData[dataIndex++];
                hexPoints[6] = cellVertData[dataIndex++];
                hexPoints[7] = cellVertData[dataIndex++];
                cellMap_[cellI] = i;
                cells_[cellI++] = cellShape(hex, hexPoints, true);
            }
            break;

            default:
                warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                dataIndex += nRead;
        }
    }

    if (debug)
    {
        Info<< "Read " << cellI << " cells;" << faceI << " faces." << endl;
    }
    cells_.setSize(cellI);
    cellMap_.setSize(cellI);
    faces_.setSize(faceI);
    faceMap_.setSize(faceI);
    lines_.setSize(lineI);
    lineMap_.setSize(lineI);
}


// Read single field and stores it on the objectRegistry.
void Foam::vtkUnstructuredReader::readField
(
    ISstream& inFile,
    objectRegistry& obj,
    const word& arrayName,
    const word& dataType,
    const label size
) const
{
    switch (vtkDataTypeNames[dataType])
    {
        case VTK_INT:
        case VTK_UINT:
        case VTK_LONG:
        case VTK_ULONG:
        case VTK_ID:
        {
            autoPtr<labelIOField> fieldVals
            (
                new labelIOField
                (
                    IOobject
                    (
                        arrayName,
                        "",
                        obj
                    ),
                    size
                )
            );
            readBlock(inFile, fieldVals().size(), fieldVals());
            regIOobject::store(fieldVals);
        }
        break;

        case VTK_FLOAT:
        case VTK_DOUBLE:
        {
            autoPtr<scalarIOField> fieldVals
            (
                new scalarIOField
                (
                    IOobject
                    (
                        arrayName,
                        "",
                        obj
                    ),
                    size
                )
            );
            readBlock(inFile, fieldVals().size(), fieldVals());
            regIOobject::store(fieldVals);
        }
        break;

        case VTK_STRING:
        {
            if (debug)
            {
                Info<< "Reading strings:" << size << endl;
            }
            autoPtr<stringIOList> fieldVals
            (
                new stringIOList
                (
                    IOobject
                    (
                        arrayName,
                        "",
                        obj
                    ),
                    size
                )
            );
            // Consume current line.
            inFile.getLine(fieldVals()[0]);
            // Read without parsing
            forAll(fieldVals(), i)
            {
                inFile.getLine(fieldVals()[i]);
            }
            regIOobject::store(fieldVals);
        }
        break;

        default:
        {
            IOWarningIn("vtkUnstructuredReader::extractCells(..)", inFile)
                << "Unhandled type " << vtkDataTypeNames[dataType] << endl
                << "Skipping " << size
                << " words." << endl;
            scalarField fieldVals;
            readBlock(inFile, size, fieldVals);
        }
        break;
    }
}


// Reads fields, stores them on the objectRegistry. Returns a list of
// read fields
Foam::wordList Foam::vtkUnstructuredReader::readFieldArray
(
    ISstream& inFile,
    objectRegistry& obj,
    const label wantedSize
) const
{
    DynamicList<word> fields;

    word dataName(inFile);
    if (debug)
    {
        Info<< "dataName:" << dataName << endl;
    }
    label numArrays(readLabel(inFile));
    if (debug)
    {
        Pout<< "numArrays:" << numArrays << endl;
    }
    for (label i = 0; i < numArrays; i++)
    {
        word arrayName(inFile);
        label numComp(readLabel(inFile));
        label numTuples(readLabel(inFile));
        word dataType(inFile);

        if (debug)
        {
            Info<< "Reading field " << arrayName
                << " of " << numTuples << " tuples of rank " << numComp << endl;
        }

        if (wantedSize != -1 && numTuples != wantedSize)
        {
            FatalIOErrorIn("vtkUnstructuredReader::readFieldArray(..)", inFile)
                << "Expected " << wantedSize << " tuples but only have "
                << numTuples << exit(FatalIOError);
        }

        readField
        (
            inFile,
            obj,
            arrayName,
            dataType,
            numTuples*numComp
        );
        fields.append(arrayName);
    }
    return fields.shrink();
}


Foam::objectRegistry& Foam::vtkUnstructuredReader::selectRegistry
(
    const parseMode readMode
)
{
    if (readMode == CELL_DATA)
    {
        return cellData_;
    }
    else if (readMode == POINT_DATA)
    {
        return pointData_;
    }
    else
    {
        return otherData_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkUnstructuredReader::vtkUnstructuredReader
(
    const objectRegistry& obr,
    ISstream& inFile
)
:
    cellData_(IOobject("cellData", obr)),
    pointData_(IOobject("pointData", obr)),
    otherData_(IOobject("otherData", obr))
{
    read(inFile);
}


void Foam::vtkUnstructuredReader::read(ISstream& inFile)
{
    inFile.getLine(header_);
    if (debug)
    {
        Info<< "Header   : " << header_ << endl;
    }
    inFile.getLine(title_);
    if (debug)
    {
        Info<< "Title    : " << title_ << endl;
    }
    inFile.getLine(dataType_);
    if (debug)
    {
        Info<< "dataType : " << dataType_ << endl;
    }

    if (dataType_ == "BINARY")
    {
        FatalIOErrorIn("vtkUnstructuredReader::read(ISstream&)", inFile)
            << "Binary reading not supported " << exit(FatalIOError);
    }

    parseMode readMode = NOMODE;
    label wantedSize = -1;


    // Temporary storage for vertices of cells.
    labelList cellVerts;

    while (inFile.good())
    {
        word tag(inFile);

        if (!inFile.good())
        {
            break;
        }

        if (debug)
        {
            Info<< "line:" << inFile.lineNumber()
                << " tag:" << tag << endl;
        }

        if (tag == "DATASET")
        {
            word geomType(inFile);
            if (debug)
            {
                Info<< "geomType : " << geomType << endl;
            }
            readMode = parseModeNames[geomType];
            wantedSize = -1;
        }
        else if (tag == "POINTS")
        {
            label nPoints(readLabel(inFile));
            points_.setSize(nPoints);    ///3);
            if (debug)
            {
                Info<< "Reading " << nPoints << " numbers representing "
                    << points_.size() << " coordinates." << endl;
            }

            word primitiveTag(inFile);
            if (primitiveTag != "float" && primitiveTag != "double")
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Expected 'float' entry but found "
                    << primitiveTag
                    << exit(FatalIOError);
            }
            forAll(points_, i)
            {
                inFile >> points_[i].x() >> points_[i].y() >> points_[i].z();
            }
        }
        else if (tag == "CELLS")
        {
            label nCells(readLabel(inFile));
            label nNumbers(readLabel(inFile));
            if (debug)
            {
                Info<< "Reading " << nCells << " cells or faces." << endl;
            }
            readBlock(inFile, nNumbers, cellVerts);
        }
        else if (tag == "CELL_TYPES")
        {
            label nCellTypes(readLabel(inFile));

            labelList cellTypes;
            readBlock(inFile, nCellTypes, cellTypes);

            if (cellTypes.size() > 0 && cellVerts.size() == 0)
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Found " << cellTypes.size()
                    << " cellTypes but no cells."
                    << exit(FatalIOError);
            }

            extractCells(inFile, cellTypes, cellVerts);
            cellVerts.clear();
        }
        else if (tag == "LINES")
        {
            label nLines(readLabel(inFile));
            label nNumbers(readLabel(inFile));
            if (debug)
            {
                Info<< "Reading " << nLines << " lines." << endl;
            }
            labelList lineVerts;
            readBlock(inFile, nNumbers, lineVerts);

            label lineI = lines_.size();
            lines_.setSize(lineI+nLines);
            lineMap_.setSize(lines_.size());

            label elemI = 0;
            for (label i = 0; i < nLines; i++)
            {
                lineMap_[lineI] = lineI;
                labelList& f = lines_[lineI];
                f.setSize(lineVerts[elemI++]);
                forAll(f, fp)
                {
                    f[fp] = lineVerts[elemI++];
                }
                lineI++;
            }
        }
        else if (tag == "POLYGONS")
        {
            // If in polydata mode

            label nFaces(readLabel(inFile));
            label nNumbers(readLabel(inFile));
            if (debug)
            {
                Info<< "Reading " << nFaces << " faces." << endl;
            }
            labelList faceVerts;
            readBlock(inFile, nNumbers, faceVerts);

            label faceI = faces_.size();
            faces_.setSize(faceI+nFaces);
            faceMap_.setSize(faces_.size());

            label elemI = 0;
            for (label i = 0; i < nFaces; i++)
            {
                faceMap_[faceI] = faceI;
                face& f = faces_[faceI];
                f.setSize(faceVerts[elemI++]);
                forAll(f, fp)
                {
                    f[fp] = faceVerts[elemI++];
                }
                faceI++;
            }
        }
        else if (tag == "POINT_DATA")
        {
            // 'POINT_DATA 24'
            readMode = POINT_DATA;
            wantedSize = points_.size();

            label nPoints(readLabel(inFile));
            if (nPoints != wantedSize)
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Reading POINT_DATA : expected " << wantedSize
                    << " but read " << nPoints << exit(FatalIOError);
            }
        }
        else if (tag == "CELL_DATA")
        {
            readMode = CELL_DATA;
            wantedSize = cells_.size()+faces_.size()+lines_.size();

            label nCells(readLabel(inFile));
            if (nCells != wantedSize)
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Reading CELL_DATA : expected "
                    << wantedSize
                    << " but read " << nCells << exit(FatalIOError);
            }
        }
        else if (tag == "FIELD")
        {
            // wantedSize already set according to type we expected to read.
            readFieldArray(inFile, selectRegistry(readMode), wantedSize);
        }
        else if (tag == "SCALARS")
        {
            string line;
            inFile.getLine(line);
            IStringStream is(line);
            word dataName(is);
            word dataType(is);
            //label numComp(readLabel(inFile));

            if (debug)
            {
                Info<< "Reading scalar " << dataName
                    << " of type " << dataType
                    << " from lookup table" << endl;
            }

            word lookupTableTag(inFile);
            if (lookupTableTag != "LOOKUP_TABLE")
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Expected tag LOOKUP_TABLE but read "
                    << lookupTableTag
                    << exit(FatalIOError);
            }

            word lookupTableName(inFile);

            readField
            (
                inFile,
                selectRegistry(readMode),
                dataName,
                dataType,
                wantedSize//*numComp
            );
        }
        else if (tag == "VECTORS" || tag == "NORMALS")
        {
            // 'NORMALS Normals float'
            string line;
            inFile.getLine(line);
            IStringStream is(line);
            word dataName(is);
            word dataType(is);
            if (debug)
            {
                Info<< "Reading vector " << dataName
                    << " of type " << dataType << endl;
            }

            objectRegistry& reg = selectRegistry(readMode);

            readField
            (
                inFile,
                reg,
                dataName,
                dataType,
                3*wantedSize
            );

            if
            (
                vtkDataTypeNames[dataType] == VTK_FLOAT
             || vtkDataTypeNames[dataType] == VTK_DOUBLE
            )
            {
                objectRegistry::iterator iter = reg.find(dataName);
                scalarField s(*dynamic_cast<const scalarField*>(iter()));
                reg.erase(iter);
                autoPtr<vectorIOField> fieldVals
                (
                    new vectorIOField
                    (
                        IOobject
                        (
                            dataName,
                            "",
                            reg
                        ),
                        s.size()/3
                    )
                );

                label elemI = 0;
                forAll(fieldVals(), i)
                {
                    fieldVals()[i].x() = s[elemI++];
                    fieldVals()[i].y() = s[elemI++];
                    fieldVals()[i].z() = s[elemI++];
                }
                regIOobject::store(fieldVals);
            }
        }
        else if (tag == "TEXTURE_COORDINATES")
        {
            // 'TEXTURE_COORDINATES TCoords 2 float'
            string line;
            inFile.getLine(line);
            IStringStream is(line);
            word dataName(is);          //"Tcoords"
            label dim(readLabel(is));
            word dataType(is);

            if (debug)
            {
                Info<< "Reading texture coords " << dataName
                    << " dimension " << dim
                    << " of type " << dataType << endl;
            }

            scalarField coords(dim*points_.size());
            readBlock(inFile, coords.size(), coords);
        }
        else if (tag == "TRIANGLE_STRIPS")
        {
            label nStrips(readLabel(inFile));
            label nNumbers(readLabel(inFile));
            if (debug)
            {
                Info<< "Reading " << nStrips << " triangle strips." << endl;
            }
            labelList faceVerts;
            readBlock(inFile, nNumbers, faceVerts);

            // Count number of triangles
            label elemI = 0;
            label nTris = 0;
            for (label i = 0; i < nStrips; i++)
            {
                label nVerts = faceVerts[elemI++];
                nTris += nVerts-2;
                elemI += nVerts;
            }


            // Store
            label faceI = faces_.size();
            faces_.setSize(faceI+nTris);
            faceMap_.setSize(faces_.size());
            elemI = 0;
            for (label i = 0; i < nStrips; i++)
            {
                label nVerts = faceVerts[elemI++];
                label nTris = nVerts-2;

                // Read first triangle
                faceMap_[faceI] = faceI;
                face& f = faces_[faceI++];
                f.setSize(3);
                f[0] = faceVerts[elemI++];
                f[1] = faceVerts[elemI++];
                f[2] = faceVerts[elemI++];
                for (label triI = 1; triI < nTris; triI++)
                {
                    faceMap_[faceI] = faceI;
                    face& f = faces_[faceI++];
                    f.setSize(3);
                    f[0] = faceVerts[elemI-1];
                    f[1] = faceVerts[elemI-2];
                    f[2] = faceVerts[elemI++];
                }
            }
        }
        else
        {
            FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                << "Unsupported tag "
                << tag << exit(FatalIOError);
        }
    }

    if (debug)
    {
        Info<< "Read points:" << points_.size()
            << " cellShapes:" << cells_.size()
            << " faces:" << faces_.size()
            << " lines:" << lines_.size()
            << nl << endl;

        Info<< "Cell fields:" << endl;
        printFieldStats<vectorIOField>(cellData_);
        printFieldStats<scalarIOField>(cellData_);
        printFieldStats<labelIOField>(cellData_);
        printFieldStats<stringIOList>(cellData_);
        Info<< nl << endl;

        Info<< "Point fields:" << endl;
        printFieldStats<vectorIOField>(pointData_);
        printFieldStats<scalarIOField>(pointData_);
        printFieldStats<labelIOField>(pointData_);
        printFieldStats<stringIOList>(pointData_);
        Info<< nl << endl;

        Info<< "Other fields:" << endl;
        printFieldStats<vectorIOField>(otherData_);
        printFieldStats<scalarIOField>(otherData_);
        printFieldStats<labelIOField>(otherData_);
        printFieldStats<stringIOList>(otherData_);
    }
}


// ************************************************************************* //

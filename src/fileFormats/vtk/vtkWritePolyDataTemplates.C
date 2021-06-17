/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "vtkWritePolyData.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DataType>
void Foam::vtkWritePolyData::writeFieldTypeValues
(
    std::ostream& os,
    const bool binary,
    const wordList& fieldNames,
    const boolList& fieldIsPointValues,
    const UPtrList<const Field<Type>>& fieldTypeValues,
    const bool writePointValues
)
{
    forAll(fieldNames, fieldi)
    {
        if
        (
            fieldIsPointValues[fieldi] == writePointValues
         && fieldTypeValues.set(fieldi)
        )
        {
            const label nCmpt = pTraits<Type>::nComponents;

            os  << fieldNames[fieldi] << ' ' << pTraits<Type>::nComponents
                << ' ' << fieldTypeValues[fieldi].size() << ' '
                << (std::is_integral<DataType>::value ? "int" : "float") << nl;

            List<DataType> data(nCmpt*fieldTypeValues[fieldi].size());
            label i = 0;
            forAll(fieldTypeValues[fieldi], fieldValuei)
            {
                for (direction cmpt = 0; cmpt < nCmpt; ++ cmpt)
                {
                    data[i ++] =
                        component
                        (
                            fieldTypeValues[fieldi][fieldValuei],
                            cmpt
                        );
                }
            }

            vtkWriteOps::write(os, binary, data);
        }
    }
}


template<class PointField, class VertexList, class LineList, class FaceList>
void Foam::vtkWritePolyData::write
(
    const fileName& file,
    const word& title,
    const bool binary,
    const PointField& points,
    const VertexList& vertices,
    const LineList& lines,
    const FaceList& faces,
    const wordList& fieldNames,
    const boolList& fieldIsPointValues,
    const UPtrList<const Field<label>>& fieldLabelValues
    #define FieldTypeValuesConstArg(Type, nullArg) \
        , const UPtrList<const Field<Type>>& field##Type##Values
    FOR_ALL_FIELD_TYPES(FieldTypeValuesConstArg)
    #undef FieldTypeValuesConstArg
)
{
    // Open the file
    std::ofstream os(file, std::ios::binary);

    // Write the header
    vtkWriteOps::writeHeader(os, binary, title);
    os << "DATASET POLYDATA" << nl;

    // Write the points
    {
        os  << "POINTS " << points.size() << " float" << nl;
        List<floatScalar> coordinates(points.size()*3);
        forAll(points, pointi)
        {
            const point& p = points[pointi];
            forAll(p, i)
            {
                coordinates[3*pointi + i] = float(p[i]);
            }
        }
        vtkWriteOps::write(os, binary, coordinates);
    }

    // Write the vertices
    if (vertices.size())
    {
        os  << "VERTICES " << vertices.size() << ' '
            << 2*vertices.size() << nl;
        labelList data(2*vertices.size());
        forAll(vertices, vertexi)
        {
            data[2*vertexi] = 1;
            data[2*vertexi + 1] = vertices[vertexi];
        }
        vtkWriteOps::write(os, binary, data);
    }

    // Write the lines
    if (lines.size())
    {
        label nLineNodes = 0;
        forAll(lines, facei)
        {
            nLineNodes += lines[facei].size();
        }
        os  << "LINES " << lines.size() << ' '
            << lines.size() + nLineNodes << nl;
        labelList data(lines.size() + nLineNodes);
        label i = 0;
        forAll(lines, linei)
        {
            data[i ++] = lines[linei].size();
            forAll(lines[linei], linePointi)
            {
                data[i ++] = lines[linei][linePointi];
            }
        }
        vtkWriteOps::write(os, binary, data);
    }

    // Write the faces
    if (faces.size())
    {
        label nFaceNodes = 0;
        forAll(faces, facei)
        {
            nFaceNodes += faces[facei].size();
        }
        os  << "POLYGONS " << faces.size() << ' '
            << faces.size() + nFaceNodes << nl;
        labelList data(faces.size() + nFaceNodes);
        label i = 0;
        forAll(faces, facei)
        {
            data[i ++] = faces[facei].size();
            forAll(faces[facei], facePointi)
            {
                data[i ++] = faces[facei][facePointi];
            }
        }
        vtkWriteOps::write(os, binary, data);
    }

    // Write the fields
    const label nPointFields = count(fieldIsPointValues, true);
    const label nFaceFields = count(fieldIsPointValues, false);
    if (nPointFields > 0)
    {
        os  << "POINT_DATA " << points.size() << nl
            << "FIELD attributes " << nPointFields << nl;
        writeFieldTypeValues<label, label>
        (
            os,
            binary,
            fieldNames,
            fieldIsPointValues,
            fieldLabelValues,
            true
        );
        #define WriteFieldTypeValues(Type, nullArg)                            \
            writeFieldTypeValues<Type, floatScalar>                            \
            (                                                                  \
                os,                                                            \
                binary,                                                        \
                fieldNames,                                                    \
                fieldIsPointValues,                                            \
                field##Type##Values,                                           \
                true                                                           \
            );
        FOR_ALL_FIELD_TYPES(WriteFieldTypeValues)
        #undef WriteFieldTypeValues
    }
    if (nFaceFields > 0)
    {
        os  << "CELL_DATA "
            << vertices.size() + lines.size() + faces.size() << nl
            << "FIELD attributes " << nFaceFields << nl;
        writeFieldTypeValues<label, label>
        (
            os,
            binary,
            fieldNames,
            fieldIsPointValues,
            fieldLabelValues,
            false
        );
        #define WriteFieldTypeValues(Type, nullArg)                            \
            writeFieldTypeValues<Type, floatScalar>                            \
            (                                                                  \
                os,                                                            \
                binary,                                                        \
                fieldNames,                                                    \
                fieldIsPointValues,                                            \
                field##Type##Values,                                           \
                false                                                          \
            );
        FOR_ALL_FIELD_TYPES(WriteFieldTypeValues)
        #undef WriteFieldTypeValues
    }
}


template
<
    class PointField,
    class VertexList,
    class LineList,
    class FaceList,
    class ... Args
>
inline void Foam::vtkWritePolyData::write
(
    const fileName& file,
    const word& title,
    const bool binary,
    const PointField& points,
    const VertexList& vertices,
    const LineList& lines,
    const FaceList& faces,
    const Args& ... args
)
{
    const label nFields = sizeof...(Args)/3;

    wordList fieldNames(nFields);
    boolList fieldIsPointValues(nFields);
    UPtrList<const Field<label>> fieldLabelValues(nFields);
    #define DeclareFieldTypeValues(Type, nullArg) \
        UPtrList<const Field<Type>> field##Type##Values(nFields);
    FOR_ALL_FIELD_TYPES(DeclareFieldTypeValues);
    #undef DeclareFieldTypeValues

    unpackFieldTypeValues
    (
        fieldNames,
        fieldIsPointValues,
        fieldLabelValues
        #define FieldTypeValuesParameter(Type, nullArg) , field##Type##Values
        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter),
        #undef FieldTypeValuesParameter
        args ...
    );

    write
    (
        file,
        title,
        binary,
        points,
        vertices,
        lines,
        faces,
        fieldNames,
        fieldIsPointValues,
        fieldLabelValues
        #define FieldTypeValuesParameter(Type, nullArg) , field##Type##Values
        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
        #undef FieldTypeValuesParameter
    );
}


// ************************************************************************* //

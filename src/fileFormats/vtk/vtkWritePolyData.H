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

Class
    Foam::vtkWritePolyData

Description
    General write functions for vtk polygonal data files

SourceFiles
    vtkWritePolyData.C

\*---------------------------------------------------------------------------*/

#ifndef vtkWritePolyData_H
#define vtkWritePolyData_H

#include "vtkWriteOps.H"
#include "fileName.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace vtkWritePolyData Declaration
\*---------------------------------------------------------------------------*/

namespace vtkWritePolyData
{

//- Helper for templated write
template<class Types, class Type>
inline void setFieldTypeValue
(
    UPtrList<const Field<Types>>& fieldTypeValues,
    const label fieldi,
    const Field<Type>& fieldTypeValue
)
{}


//- Helper for templated write
template<class Type>
inline void setFieldTypeValue
(
    UPtrList<const Field<Type>>& fieldTypeValues,
    const label fieldi,
    const Field<Type>& fieldTypeValue
)
{
    fieldTypeValues.set(fieldi, &fieldTypeValue);
}


//- Helper for templated write
inline void unpackFieldTypeValues
(
    wordList& fieldNames,
    boolList& fieldIsPointValues,
    UPtrList<const Field<label>>& fieldLabelValues
    #define FieldTypeValuesNonConstArg(Type, nullArg) \
    , UPtrList<const Field<Type>>& field##Type##Values
    FOR_ALL_FIELD_TYPES(FieldTypeValuesNonConstArg)
    #undef FieldTypeValuesNonConstArg
)
{}


//- Helper for templated write
template<class Type, class ... Args>
inline void unpackFieldTypeValues
(
    wordList& fieldNames,
    boolList& fieldIsPointValues,
    UPtrList<const Field<label>>& fieldLabelValues
    #define FieldTypeValuesNonConstArg(Type, nullArg) \
    , UPtrList<const Field<Type>>& field##Type##Values
    FOR_ALL_FIELD_TYPES(FieldTypeValuesNonConstArg),
    #undef FieldTypeValuesNonConstArg
    const word& fieldName,
    const bool fieldIsPointValue,
    const Field<Type>& fieldTypeValue,
    Args& ... args
)
{
    const label fieldi = fieldNames.size() - 1 - sizeof...(Args)/3;

    fieldNames[fieldi] = fieldName;
    fieldIsPointValues[fieldi] = fieldIsPointValue;
    setFieldTypeValue(fieldLabelValues, fieldi, fieldTypeValue);
    #define SetFieldTypeValue(Type, nullArg) \
        setFieldTypeValue(field##Type##Values, fieldi, fieldTypeValue);
    FOR_ALL_FIELD_TYPES(SetFieldTypeValue);
    #undef SetFieldTypeValue

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
}


//- Write the field values out for a type
template<class Type, class DataType>
void writeFieldTypeValues
(
    std::ostream& os,
    const bool binary,
    const wordList& fieldNames,
    const boolList& fieldIsPointValues,
    const UPtrList<const Field<Type>>& fieldTypeValues,
    const bool writePointValues
);


//- Write VTK polygonal data to a file. Takes a PtrList of fields of labels and
//  of every primitive type. Each PtrList should be the length of the total
//  number of fields and only one PtrList should be non-null for each field
//  index.
template<class PointField, class VertexList, class LineList, class FaceList>
void write
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
);


//- Write VTK polygonal data to a file. Takes any number of name,
//  isPointValues, values arguments at the end. E.g.;
//
//  write
//  (
//      // Output options
//      "myPolyData.vtk", "myPolyData", false,
//
//      // Geometry
//      pp.localPoints(), labelList(), labelListList(), pp.localFaces(),
//
//      // Fields
//      "faceIDs", true, pp.addressing(),
//      "facePressures", false, Field<scalar>(pp.size(), ...),
//      "pointVelocities", true, Field<vector>(pp.nPoints(), ...)
//  );
//
template
<
    class PointField,
    class VertexList,
    class LineList,
    class FaceList,
    class ... Args
>
inline void write
(
    const fileName& file,
    const word& title,
    const bool binary,
    const PointField& points,
    const VertexList& vertices,
    const LineList& lines,
    const FaceList& faces,
    const Args& ... args
);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtkWritePolyDataTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

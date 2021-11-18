/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "rawSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rawSurfaceWriter, 0);
    addToRunTimeSelectionTable(surfaceWriter, rawSurfaceWriter, wordDict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::rawSurfaceWriter::writeLocation
(
    Ostream& os,
    const pointField& points,
    const label pointi
)
{
    const point& pt = points[pointi];
    os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << ' ';
}


inline void Foam::rawSurfaceWriter::writeLocation
(
    Ostream& os,
    const pointField& points,
    const faceList& faces,
    const label facei
)
{
    const point& ct = faces[facei].centre(points);
    os  << ct.x() << ' ' << ct.y() << ' ' << ct.z() << ' ';
}


template<class Type>
void Foam::rawSurfaceWriter::writeHeader
(
    Ostream& os,
    const word& fieldName
)
{
    const label nCmpt = pTraits<Type>::nComponents;

    for (direction cmpt = 0; cmpt < nCmpt; ++ cmpt)
    {
        const bool separator =
            !fieldName.empty()
         && strlen(pTraits<Type>::componentNames[cmpt]) > 0;

        if (cmpt) os  << token::SPACE;
        os  << fieldName << (separator ? "_" : "")
            << pTraits<Type>::componentNames[cmpt];
    }
}


template<class Type>
inline void Foam::rawSurfaceWriter::writeData
(
    Ostream& os,
    const Type& v
)
{
    const label nCmpt = pTraits<Type>::nComponents;

    for (direction cmpt = 0; cmpt < nCmpt; ++ cmpt)
    {
        if (cmpt) os  << token::SPACE;
        os  << component(v, cmpt);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rawSurfaceWriter::rawSurfaceWriter
(
    const IOstream::streamFormat writeFormat
)
:
    surfaceWriter(writeFormat),
    writeCompression_(IOstream::UNCOMPRESSED)
{}


Foam::rawSurfaceWriter::rawSurfaceWriter(const dictionary& optDict)
:
    surfaceWriter(optDict),
    writeCompression_(IOstream::UNCOMPRESSED)
{
    if (optDict.found("compression"))
    {
        writeCompression_ =
            IOstream::compressionEnum(optDict.lookup("compression"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rawSurfaceWriter::~rawSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rawSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const wordList& fieldNames,
    const bool writePointValues
    #define FieldTypeValuesConstArg(Type, nullArg) \
        , const UPtrList<const Field<Type>>& field##Type##Values
    FOR_ALL_FIELD_TYPES(FieldTypeValuesConstArg)
    #undef FieldTypeValuesConstArg
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os
    (
        outputDir/surfaceName + ".raw",
        IOstream::ASCII,
        IOstream::currentVersion,
        writeCompression_
    );

    // Get the number of values
    label nValues = 0;
    #define GetNValues(Type, nullArg)                \
        if (field##Type##Values.set(0))              \
        {                                            \
            nValues = field##Type##Values[0].size(); \
        }
    FOR_ALL_FIELD_TYPES(GetNValues);
    #undef GetNValues

    // Top header
    os << "# " << (writePointValues ? "POINT_DATA " : "FACE_DATA ")
        << nValues << nl;

    // Column headers
    os << "# ";
    writeHeader<vector>(os, word::null);
    forAll(fieldNames, fieldi)
    {
        os << token::SPACE;

        #define WriteTypeHeader(Type, nullArg)             \
            if (field##Type##Values.set(fieldi))           \
            {                                              \
                writeHeader<Type>(os, fieldNames[fieldi]); \
            }
        FOR_ALL_FIELD_TYPES(WriteTypeHeader);
        #undef WriteTypeHeader
    }
    os << nl;

    // Write the values
    #define WriteTypeValues(Type, nullArg)                 \
        if (field##Type##Values.set(fieldi))               \
        {                                                  \
            writeData(os, field##Type##Values[fieldi][i]); \
        }
    if (writePointValues)
    {
        for (label i = 0; i < nValues; ++ i)
        {
            writeLocation(os, points, i);
            forAll(fieldNames, fieldi)
            {
                if (fieldi) os << token::SPACE;
                FOR_ALL_FIELD_TYPES(WriteTypeValues);
            }
            os << nl;
        }
    }
    else
    {
        for (label i = 0; i < nValues; ++ i)
        {
            writeLocation(os, points, faces, i);
            forAll(fieldNames, fieldi)
            {
                if (fieldi) os << token::SPACE;
                FOR_ALL_FIELD_TYPES(WriteTypeValues);
            }
            os << nl;
        }
    }
    #undef WriteTypeValues
}


// ************************************************************************* //

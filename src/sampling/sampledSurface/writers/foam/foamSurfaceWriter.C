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

#include "foamSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(foamSurfaceWriter, 0);
    addToRunTimeSelectionTable(surfaceWriter, foamSurfaceWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, foamSurfaceWriter, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamSurfaceWriter::~foamSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamSurfaceWriter::write
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
    const fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (debug)
    {
        Info<< "Writing geometry to " << surfaceDir << endl;
    }

    auto stream = [&](const fileName& file)
    {
        return
            OFstream
            (
                file,
                writeFormat_,
                IOstream::currentVersion,
                writeCompression_
            );
    };

    // Points
    stream(surfaceDir/"points")() << points;

    // Faces
    stream(surfaceDir/"faces")() << faces;

    // Face centers. Not really necessary but very handy when reusing as inputs
    // for e.g. timeVaryingMapped bc.
    pointField faceCentres(faces.size(), Zero);
    forAll(faces, facei)
    {
        faceCentres[facei] = faces[facei].centre(points);
    }
    stream(surfaceDir/"faceCentres")() << faceCentres;

    // Fields
    forAll(fieldNames, fieldi)
    {
        if (debug)
        {
            Info<< "Writing field " << fieldNames[fieldi] << " to "
                << surfaceDir << endl;
        }

        #define WriteFieldType(Type, nullArg)          \
            if (field##Type##Values.set(fieldi))       \
            {                                          \
                const fileName valuesDir               \
                (                                      \
                    surfaceDir                         \
                   /(                                  \
                       word(pTraits<Type>::typeName)   \
                     + word(Field<Type>::typeName)     \
                    )                                  \
                );                                     \
                                                       \
                if (!isDir(valuesDir))                 \
                {                                      \
                    mkDir(valuesDir);                  \
                }                                      \
                                                       \
                stream(valuesDir/fieldNames[fieldi])() \
                    << field##Type##Values[fieldi];    \
            }
        FOR_ALL_FIELD_TYPES(WriteFieldType);
        #undef WriteFieldType
    }
}


// ************************************************************************* //

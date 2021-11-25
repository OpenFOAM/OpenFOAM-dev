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

#include "ensightSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "ensightPartFaces.H"
#include "ensightPTraits.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightSurfaceWriter, 0);
    addToRunTimeSelectionTable(surfaceWriter, ensightSurfaceWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, ensightSurfaceWriter, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightSurfaceWriter::~ensightSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightSurfaceWriter::write
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

    // const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;

    OFstream osCase(outputDir/surfaceName + ".case");
    ensightGeoFile osGeom
    (
        outputDir/surfaceName + ".000.mesh",
        writeFormat_
    );

    if (debug)
    {
        Info<< "Writing case file to " << osCase.name() << endl;
    }

    osCase
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl;

    osCase
        << "GEOMETRY" << nl
        << "model:        1     " << osGeom.name().name() << nl
        << nl;

    osCase
        << "VARIABLE" << nl;
    forAll(fieldNames, fieldi)
    {
        #define WriteTypeCase(Type, nullArg)                            \
            if (field##Type##Values.set(fieldi))                        \
            {                                                           \
                osCase                                                  \
                    << ensightPTraits<Type>::typeName << " per "        \
                    << word(writePointValues ? "node:" : "element:")    \
                    << setw(10) << 1 << "       " << fieldNames[fieldi] \
                    << "       " << surfaceName.c_str() << ".***."      \
                    << fieldNames[fieldi] << nl;                        \
            }
        FOR_ALL_FIELD_TYPES(WriteTypeCase);
        #undef WriteTypeCase
    }
    osCase
        << nl;

    osCase
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << timeValue << nl
        << nl;

    ensightPartFaces ensPart(0, osGeom.name().name(), points, faces, true);
    osGeom << ensPart;

    forAll(fieldNames, fieldi)
    {
        #define WriteTypeValues(Type, nullArg)                        \
            if (field##Type##Values.set(fieldi))                      \
            {                                                         \
                ensightFile osField                                   \
                (                                                     \
                    outputDir/surfaceName                             \
                  + ".000."                                           \
                  + fieldNames[fieldi],                               \
                    writeFormat_                                      \
                );                                                    \
                osField.writeKeyword(ensightPTraits<Type>::typeName); \
                ensPart.writeField                                    \
                (                                                     \
                    osField,                                          \
                    field##Type##Values[fieldi],                      \
                    writePointValues                                  \
                );                                                    \
            }
        FOR_ALL_FIELD_TYPES(WriteTypeValues);
        #undef WriteTypeValues
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::surfaceWriter

Description
    Base class for surface writers

SourceFiles
    surfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWriter_H
#define surfaceWriter_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "pointField.H"
#include "faceList.H"
#include "fileName.H"
#include "setWriter.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceWriter
{
protected:

    // Protected Data

        //- Write format
        IOstream::streamFormat writeFormat_;

        //- Write compression
        IOstream::compressionType writeCompression_;


public:

    //- Runtime type information
    TypeName("surfaceWriter");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            word,
            (
                const IOstream::streamFormat writeFormat,
                const IOstream::compressionType writeCompression
            ),
            (writeFormat, writeCompression)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            dict,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Select given write options
        static autoPtr<surfaceWriter> New
        (
            const word& writeType,
            const IOstream::streamFormat writeFormat,
            const IOstream::compressionType writeCompression
        );

        //- Select given a dictionary
        static autoPtr<surfaceWriter> New
        (
            const word& writeType,
            const dictionary& dict
        );


    // Constructors

        //- Construct given write options
        surfaceWriter
        (
            const IOstream::streamFormat writeFormat,
            const IOstream::compressionType writeCompression
        );

        //- Construct from dictionary
        surfaceWriter(const dictionary& dict);


    //- Destructor
    virtual ~surfaceWriter();


    // Member Functions

        //- Write fields for a single surface to file.
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const wordList& fieldNames,     // names of fields
            const bool writePointValues
            #define FieldTypeValuesConstArg(Type, nullArg) \
                , const UPtrList<const Field<Type>>& field##Type##Values
            FOR_ALL_FIELD_TYPES(FieldTypeValuesConstArg)
            #undef FieldTypeValuesConstArg
        ) const = 0;

        //- Write fields for a single surface to file.
        virtual void write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const wordList& fieldNames,     // names of fields
            const bool writePointValues
            #define FieldTypeValuesConstArg(Type, nullArg) \
                , const PtrList<Field<Type>>& field##Type##Values
            FOR_ALL_FIELD_TYPES(FieldTypeValuesConstArg)
            #undef FieldTypeValuesConstArg
        ) const
        {
            write
            (
                outputDir,
                surfaceName,
                points,
                faces,
                fieldNames,
                writePointValues
                #define FieldTypeValuesParameter(Type, nullArg)            \
                    , reinterpret_cast<const UPtrList<const Field<Type>>&> \
                    (field##Type##Values)
                FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                #undef FieldTypeValuesParameter
            );
        }

        //- Write the surface geometry only.
        void inline write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces
        ) const
        {
            write
            (
                outputDir,
                surfaceName,
                points,
                faces,
                wordList(),
                false
                #define FieldTypeValuesParameter(Type, nullArg) \
                    , UPtrList<const Field<Type>>()
                FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                #undef FieldTypeValuesParameter
            );
        }

        //- Write fields for a single surface to file. For use in code where
        //  the fields that are to be written are known. Takes any number of
        //  name, values arguments at the end. E.g.:
        //
        //  write
        //  (
        //      // Output options
        //      "myDirectory", "mySurface",
        //
        //      // Geometry
        //      pp.localPoints(), pp.localFaces(),
        //
        //      // Fields
        //      "p", Field<scalar>(pp.size(), ...),
        //      "U", Field<vector>(pp.size(), ...)
        //  );
        //
        template<class ... Args>
        void inline write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const bool writePointValues,
            const Args& ... args
        ) const
        {
            wordList fieldNames;
            #define DeclareFieldTypeValues(Type, nullArg) \
                UPtrList<const Field<Type>> field##Type##Values;
            FOR_ALL_FIELD_TYPES(DeclareFieldTypeValues);
            #undef DeclareFieldTypeValues

            setWriter::unpackTypeValueSets
            (
                fieldNames
                #define FieldTypeValuesParameter(Type, nullArg) \
                    , field##Type##Values
                FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter),
                #undef FieldTypeValuesParameter
                args ...
            );

            write
            (
                outputDir,
                surfaceName,
                points,
                faces,
                fieldNames,
                writePointValues
                #define FieldTypeValuesParameter(Type, nullArg) \
                    , field##Type##Values
                FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                #undef FieldTypeValuesParameter
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

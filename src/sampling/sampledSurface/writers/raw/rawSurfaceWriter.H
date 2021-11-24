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

Class
    Foam::rawSurfaceWriter

Description
    A surfaceWriter for raw output.

SourceFiles
    rawSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef rawSurfaceWriter_H
#define rawSurfaceWriter_H

#include "surfaceWriter.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rawSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class rawSurfaceWriter
:
    public surfaceWriter
{
    // Private Data

        //- Output compression, defaults to uncompressed
        IOstream::compressionType writeCompression_;


    // Private Member Functions

        //- Write a point location
        static inline void writeLocation
        (
            Ostream&,
            const pointField&,
            const label pointi
        );

        //- Write a face location
        static inline void writeLocation
        (
            Ostream&,
            const pointField&,
            const faceList&,
            const label facei
        );

        //- Write the column header for a given field
        template<class Type>
        static void writeHeader
        (
            Ostream&,
            const word& fieldName
        );

        //- Write a piece of data
        template<class Type>
        static inline void writeData(Ostream&, const Type&);


public:

    //- Runtime type information
    TypeName("raw");


    // Constructors

        //- Construct given write format
        rawSurfaceWriter(const IOstream::streamFormat writeFormat);

        //- Construct with output options
        rawSurfaceWriter(const dictionary& optDict);


    //- Destructor
    virtual ~rawSurfaceWriter();


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
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

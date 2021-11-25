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
#include "rawSetWriter.H"

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
public:

    //- Runtime type information
    TypeName("raw");


    // Constructors

        //- Inherit constructors
        using surfaceWriter::surfaceWriter;


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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::ensightSetWriter

Description

SourceFiles
    ensightSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef ensightSetWriter_H
#define ensightSetWriter_H

#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightSetWriter Declaration
\*---------------------------------------------------------------------------*/

class ensightSetWriter
:
    public setWriter
{
public:

    //- Runtime type information
    TypeName("ensight");


    // Constructors

        //- Inherit constructors
        using setWriter::setWriter;

        //- Construct and return a clone
        virtual autoPtr<setWriter> clone() const
        {
            return autoPtr<setWriter>(new ensightSetWriter(*this));
        }


    //- Destructor
    virtual ~ensightSetWriter();


    // Member Functions

        //- Write a coordSet and associated data
        virtual void write
        (
            const fileName& outputDir,
            const fileName& setName,
            const coordSet& set,
            const wordList& valueSetNames
            #define TypeValueSetsConstArg(Type, nullArg) \
                , const UPtrList<const Field<Type>>& Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
            #undef TypeValueSetsConstArg
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

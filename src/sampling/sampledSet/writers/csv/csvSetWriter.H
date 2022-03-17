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
    Foam::csvSetWriter

Description
    Write set in csv format

SourceFiles
    csvSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef csvSetWriter_H
#define csvSetWriter_H

#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class csvSetWriter Declaration
\*---------------------------------------------------------------------------*/

class csvSetWriter
:
    public setWriter
{
private:

    // Private Data

        //- Separator character
        const char separator_;


protected:

    // Protected Member Functions

        //- Write a value separator
        virtual void writeValueSeparator(Ostream& os) const;

        //- Write a segment separator
        virtual void writeSegmentSeparator(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("csv");


    // Constructors

        //- Construct given write options
        csvSetWriter
        (
            const IOstream::streamFormat writeFormat,
            const IOstream::compressionType writeCompression
        );

        //- Construct from dictionary
        csvSetWriter(const dictionary& dict);


    //- Destructor
    virtual ~csvSetWriter();


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

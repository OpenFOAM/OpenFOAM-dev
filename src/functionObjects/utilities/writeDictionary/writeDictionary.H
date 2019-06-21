/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::functionObjects::writeDictionary

Description
    Writes dictionaries on start-up and on change.

SourceFiles
    writeDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeDictionary_H
#define functionObjects_writeDictionary_H

#include "functionObject.H"
#include "wordList.H"
#include "SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeDictionary Declaration
\*---------------------------------------------------------------------------*/

class writeDictionary
:
    public functionObject
{
    // Private Data

        //- Reference to the database
        const objectRegistry& obr_;

        //- Names of dictionaries to monitor
        wordList dictNames_;

        //- List of changed dictionaries (only those registered to database)
        List<SHA1Digest> digests_;


    // Private Member Functions

        //- Helper function to write the dictionary if found at location
        bool tryDirectory
        (
            const label dictI,
            const word& location,
            bool& firstDict
        );


public:

    //- Runtime type information
    TypeName("writeDictionary");


    // Constructors

        //- Construct from Time and dictionary
        writeDictionary
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        writeDictionary(const writeDictionary&) = delete;


    //- Destructor
    virtual ~writeDictionary();


    // Member Functions

        //- Read the writeDictionary data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the selected dictionaries
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeDictionary&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

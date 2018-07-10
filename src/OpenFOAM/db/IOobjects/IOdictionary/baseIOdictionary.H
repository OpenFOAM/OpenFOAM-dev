/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::baseIOdictionary

Description
    baseIOdictionary is derived from dictionary and IOobject to give the
    dictionary automatic IO functionality via the objectRegistry.
    To facilitate IO,
    IOdictionary is provided with a constructor from IOobject and writeData and
    write functions.

SourceFiles
    baseIOdictionary.C
    baseIOdictionaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef baseIOdictionary_H
#define baseIOdictionary_H

#include "dictionary.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class baseIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class baseIOdictionary
:
    public regIOobject,
    public dictionary
{
    // Private data

        static bool writeDictionaries;

public:

    TypeName("dictionary");


    // Constructors

        //- Construct given an IOobject
        baseIOdictionary(const IOobject&);

        //- Construct given an IOobject and dictionary
        baseIOdictionary(const IOobject&, const dictionary&);

        //- Construct given an IOobject and Istream
        baseIOdictionary(const IOobject&, Istream&);


    //- Destructor
    virtual ~baseIOdictionary();


    // Member functions

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const = 0;

        //- Name function is needed to disambiguate those inherited
        //  from regIOobject and dictionary
        const word& name() const;

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const = 0;


    // Member operators

        //- Assignment of other baseIOdictionary's entries to this
        //  baseIOdictionary
        void operator=(const baseIOdictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

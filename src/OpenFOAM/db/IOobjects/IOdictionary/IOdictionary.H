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
    Foam::IOdictionary

Description
    IOdictionary is derived from dictionary and IOobject to give the
    dictionary automatic IO functionality via the objectRegistry.  To facilitate
    IO, IOdictionary is provided with a constructor from IOobject and writeData
    and write functions.

SourceFiles
    IOdictionary.C
    IOdictionaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef IOdictionary_H
#define IOdictionary_H

#include "dictionary.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IOdictionary Declaration
\*---------------------------------------------------------------------------*/

class IOdictionary
:
    public regIOobject,
    public dictionary
{

protected:

    // Protected constructors

        //- Construct given an IOobject, supply wanted typeName
        IOdictionary(const IOobject& io, const word& wantedType);


public:

    TypeName("dictionary");


    // Public static data

        static bool writeDictionaries;


    // Constructors

        //- Construct given an IOobject
        IOdictionary(const IOobject&);

        //- Construct given an IOobject and dictionary
        IOdictionary(const IOobject&, const dictionary&);

        //- Construct given an IOobject and Istream
        IOdictionary(const IOobject&, Istream&);

        //- Copy constructor
        IOdictionary(const IOdictionary&);

        //- Move constructor
        IOdictionary(IOdictionary&&);


    //- Destructor
    virtual ~IOdictionary();


    // Member Functions

        using regIOobject::name;

        //- Return true as object is global, i.e. same for all processors
        virtual bool global() const
        {
            return true;
        }

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;


    // Member Operators

        //- Assignment of dictionary entries only
        void operator=(const IOdictionary&);

        //- Move assignment of dictionary entries only
        void operator=(IOdictionary&&);
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<IOdictionary>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::IOdictionary

Description
    IOdictionary is derived from dictionary and IOobject to give the dictionary
    automatic IO functionality via the objectRegistry.  To facilitate IO,
    IOdictionary is provided with a constructor from IOobject and writeData and
    write functions.

SourceFiles
    IOdictionary.C
    IOdictionaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef IOdictionary_H
#define IOdictionary_H

#include "baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOdictionary Declaration
\*---------------------------------------------------------------------------*/

class IOdictionary
:
    public baseIOdictionary
{

public:

    // Constructors

        //- Construct given an IOobject
        IOdictionary(const IOobject&);

        //- Construct given an IOobject and dictionary
        IOdictionary(const IOobject&, const dictionary&);

        //- Construct given an IOobject and Istream
        IOdictionary(const IOobject&, Istream&);


    //- Destructor
    virtual ~IOdictionary();


    // Member functions

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }
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

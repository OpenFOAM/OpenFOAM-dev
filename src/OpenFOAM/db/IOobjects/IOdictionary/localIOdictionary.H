/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::localIOdictionary

Description
    localIOdictionary is derived from IOdictionary but excludes parallel
    master reading.

SourceFiles
    localIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef localIOdictionary_H
#define localIOdictionary_H

#include "baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class localIOdictionary
:
    public baseIOdictionary
{

public:

    // Constructors

        //- Construct given an IOobject
        localIOdictionary(const IOobject& io);

        //- Construct given an IOobject and dictionary
        localIOdictionary(const IOobject& io, const dictionary& dict);

        //- Construct given an IOobject and Istream
        localIOdictionary(const IOobject& io, Istream& is);

        //- Construct given an IOobject, supply wanted typeName
        localIOdictionary(const IOobject& io, const word& wantedType);


    //- Destructor
    virtual ~localIOdictionary();


    // Member functions

        //- Is object global
        virtual bool global() const
        {
            return false;
        }

        //- Return complete path + object name if the file exists
        //  in the case otherwise null
        virtual fileName filePath() const
        {
            // Use default (local only) search strategy
            return localFilePath(type());
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

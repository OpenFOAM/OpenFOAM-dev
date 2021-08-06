/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    localIOdictionary derived from IOdictionary with global set false to
    disable parallel master reading.

    Used to read GeometricFields.

SourceFiles
    localIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef localIOdictionary_H
#define localIOdictionary_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class localIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class localIOdictionary
:
    public IOdictionary
{

public:

    // Constructors

        //- Construct given an IOobject and actual type name
        localIOdictionary(const IOobject& io);

        //- Construct given an IOobject and dictionary
        localIOdictionary(const IOobject&, const dictionary&);

        //- Construct given an IOobject and actual type name
        localIOdictionary(const IOobject& io, const word& actualType);


    // Member Functions

        //- Return false as object is not global,
        //  i.e. different for all processors
        virtual bool global() const
        {
            return false;
        }
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<localIOdictionary>()
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

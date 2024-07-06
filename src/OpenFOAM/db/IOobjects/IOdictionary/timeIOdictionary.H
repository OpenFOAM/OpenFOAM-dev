/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::timeIOdictionary

Description
    timeIOdictionary derived from IOdictionary with globalFile set false to
    enable writing to processor time directories.

    Used for time-dependent global data written to the \<time\>/uniform
    directories.

SourceFiles
    timeIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef timeIOdictionary_H
#define timeIOdictionary_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class timeIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class timeIOdictionary
:
    public IOdictionary
{

public:

    // Constructors

        //- Construct given an IOobject and actual type name
        timeIOdictionary(const IOobject& io);

        //- Construct given an IOobject and dictionary
        timeIOdictionary(const IOobject&, const dictionary&);


    // Member Functions

        //- Return false as the object is global, i.e. same for all processors
        //  but written to the processor time directory
        virtual bool globalFile() const
        {
            return false;
        }
};


//- Trait for obtaining global write status
template<>
struct typeGlobalFile<timeIOdictionary>
{
    static const bool global = false;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

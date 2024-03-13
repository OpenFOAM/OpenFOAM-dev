/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::controlIOdictionary

Description
    controlDict specific IOdictionary to provide automatic read-update for Time

SourceFiles
    controlIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef controlIOdictionary_H
#define controlIOdictionary_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;

/*---------------------------------------------------------------------------*\
                     Class controlIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class controlIOdictionary
:
    public IOdictionary
{
    // Private Data

        //- Reference to Time to allow readDict() to be called
        //  if the controlDict file changes
        Time& time_;


public:

    // Constructors

        //- Construct for Time
        controlIOdictionary
        (
            const IOobject& io,
            Time& time
        );

        //- Construct for Time
        controlIOdictionary
        (
            const IOobject& io,
            const dictionary& dict,
            Time& time
        );

        //- Disallow default bitwise copy construction
        controlIOdictionary(const controlIOdictionary&) = delete;


    // Member Functions

        // Read

            //- Read the controlDict and update Time
            bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const controlIOdictionary&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

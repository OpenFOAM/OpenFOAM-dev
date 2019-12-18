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
    Foam::demandDrivenEntry

Description
    Class for demand-driven dictionary entries

    Holds a reference to a dictionary, which is then queried if the value
    is requested and has not already been cached

SourceFiles
    demandDrivenEntry.C
    demandDrivenEntryI.H

\*---------------------------------------------------------------------------*/

#ifndef demandDrivenEntry_H
#define demandDrivenEntry_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class demandDrivenEntry Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class demandDrivenEntry
{
    // Private Data

        //- Reference to the dictionary
        const dictionary& dict_;

        //- Keyword to look up
        const word keyword_;

        //- Value
        mutable Type value_;

        //- Flag to say that the value has been stored
        mutable bool stored_;


public:

    //- Constructors

        //- Construct from dictionary and value - cannot be re-read
        demandDrivenEntry
        (
            const dictionary& dict,
            const Type& value
        );


        //- Construct from dictionary and keyword
        demandDrivenEntry
        (
            const dictionary& dict,
            const word& keyword
        );


        //- Construct from dictionary, keyword and default value
        demandDrivenEntry
        (
            const dictionary& dict,
            const word& keyword,
            const Type& defaultValue,
            const bool readIfPresent = true
        );

        //- Copy constructor
        demandDrivenEntry(const demandDrivenEntry& dde);


    // Public Member Functions

        //- Initialise
        inline void initialise() const;

        //- Return the value
        inline const Type& value() const;

        //- Set the value
        inline void setValue(const Type& value);

        //- Reset the demand-driven entry
        inline void reset();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "demandDrivenEntryI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "demandDrivenEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

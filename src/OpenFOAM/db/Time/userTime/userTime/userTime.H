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
    Foam::userTime

Description
    An abstract class for the user time description

SourceFiles
    userTime.C

\*---------------------------------------------------------------------------*/

#ifndef userTime_H
#define userTime_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace userTimes
{

/*---------------------------------------------------------------------------*\
                         Class userTime Declaration
\*---------------------------------------------------------------------------*/

class userTime
{

public:

    //- Runtime type information
    TypeName("userTime");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        userTime,
        dictionary,
        (const dictionary& controlDict),
        (controlDict)
    );


    // Constructors

        //- Construct from controlDict
        userTime(const dictionary& controlDict);


    // Selector

        static autoPtr<userTime> New(const dictionary& controlDict);


    //- Destructor
    virtual ~userTime();


    // Member Functions

        //- Return the userTime sub-dictionary
        static const dictionary& dict(const dictionary& controlDict);

        //- Convert the user-time (e.g. CA deg) to real-time (s).
        virtual scalar userTimeToTime(const scalar tau) const = 0;

        //- Convert the real-time (s) into user-time (e.g. CA deg)
        virtual scalar timeToUserTime(const scalar t) const = 0;

        //- Return user-time unit
        virtual word unitName() const = 0;

        //- Return the user-time unit conversion
        virtual const unitConversion& units() const = 0;

        //- Read the controlDict and set all the parameters
        virtual bool read(const dictionary& controlDict) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace userTimes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::userTimes::real

Description
    A simple userTime which simply returns the real-time is seconds.

SourceFiles
    realTime.C

\*---------------------------------------------------------------------------*/

#ifndef realTime_H
#define realTime_H

#include "userTime.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace userTimes
{

/*---------------------------------------------------------------------------*\
                         Class real Declaration
\*---------------------------------------------------------------------------*/

class real
:
    public userTime
{
public:

    //- Runtime type information
    TypeName("real");


    // Constructors

        //- Construct from controlDict
        real(const dictionary& controlDict);


    //- Destructor
    virtual ~real();


    // Member Functions

        //- Return tau (s)
        virtual scalar userTimeToTime(const scalar tau) const;

        //- Return t (s)
        virtual scalar timeToUserTime(const scalar t) const;

        //- Return real-time unit name (s)
        virtual word unitName() const;

        //- Return the real-time unit conversion
        virtual const unitConversion& units() const;

        //- Read the controlDict and set all the parameters
        virtual bool read(const dictionary& controlDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace userTimes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

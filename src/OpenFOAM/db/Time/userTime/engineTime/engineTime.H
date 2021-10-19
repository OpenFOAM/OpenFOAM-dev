/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::userTimes::engine

Description
    A userTime based on engine crank-angle degrees.

SourceFiles
    engineTime.C

\*---------------------------------------------------------------------------*/

#ifndef engineTime_H
#define engineTime_H

#include "userTime.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace userTimes
{

/*---------------------------------------------------------------------------*\
                         Class engine Declaration
\*---------------------------------------------------------------------------*/

class engine
:
    public userTime
{
    // Private Data

        //- The engine RPM
        dimensionedScalar rpm_;


public:

    //- Runtime type information
    TypeName("engine");


    // Constructors

        //- Construct from controlDict
        engine(const dictionary& controlDict);


    //- Destructor
    virtual ~engine();


    // Member Functions

        //- Return the RPM
        const dimensionedScalar& rpm() const;

        //- Return the theta crank-angle is s
        virtual scalar userTimeToTime(const scalar theta) const;

        //- Return the time t in crank-angle
        virtual scalar timeToUserTime(const scalar t) const;

        //- Return real-time unit (s)
        virtual word unit() const;

        //- Read the controlDict and set all the parameters
        virtual bool read(const dictionary& controlDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace userTimes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

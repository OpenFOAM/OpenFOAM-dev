/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::TimeState

Description
    The time value with time-stepping information, user-defined remapping, etc.

SourceFiles
    TimeState.C

\*---------------------------------------------------------------------------*/

#ifndef TimeState_H
#define TimeState_H

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;


/*---------------------------------------------------------------------------*\
                          Class TimeState Declaration
\*---------------------------------------------------------------------------*/

class TimeState
:
    public dimensionedScalar
{
    // Private Member Data

        label timeIndex_;
        scalar deltaT_;
        scalar deltaTSave_;
        scalar deltaT0_;
        bool deltaTchanged_;
        label writeTimeIndex_;
        bool writeTime_;


public:

    //- Declare friendship with the Time class
    friend class Time;


    // Constructors

        TimeState();


    // Member Functions

        // Access

            //- Return current time index
            inline label timeIndex() const;

            //- Return time step value
            inline scalar deltaTValue() const;

            //- Return old time step value
            inline scalar deltaT0Value() const;

            //- Return time step
            inline dimensionedScalar deltaT() const;

            //- Return old time step
            inline dimensionedScalar deltaT0() const;


        // Check

            //- Return true if this is a write time
            inline bool writeTime() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TimeStateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

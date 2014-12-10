/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

/*---------------------------------------------------------------------------*\
                           Class TimeState Declaration
\*---------------------------------------------------------------------------*/


class TimeState
:
    public dimensionedScalar
{

protected:

        label timeIndex_;
        scalar deltaT_;
        scalar deltaTSave_;
        scalar deltaT0_;
        bool deltaTchanged_;

        label outputTimeIndex_;
        //- Is outputTime because of primary?
        bool primaryOutputTime_;

        label secondaryOutputTimeIndex_;
        //- Is outputTime because of secondary?
        bool secondaryOutputTime_;

        bool outputTime_;

public:

    // Constructors

        TimeState();


    //- Destructor
    virtual ~TimeState();


    // Member functions

        // Access

            //- Convert the user-time (e.g. CA deg) to real-time (s).
            virtual scalar userTimeToTime(const scalar theta) const;

            //- Convert the real-time (s) into user-time (e.g. CA deg)
            virtual scalar timeToUserTime(const scalar t) const;

            //- Return current time value
            scalar timeOutputValue() const;

            //- Return current time index
            label timeIndex() const;

            //- Return time step value
            inline scalar deltaTValue() const
            {
                return deltaT_;
            }

            //- Return old time step value
            inline scalar deltaT0Value() const
            {
                return deltaT0_;
            }

            //- Return time step
            dimensionedScalar deltaT() const;

            //- Return old time step
            dimensionedScalar deltaT0() const;


        // Check

            //- Return true if this is an output time (primary or secondary)
            bool outputTime() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

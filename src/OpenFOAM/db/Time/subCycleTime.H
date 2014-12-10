/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::subCycleTime

Description
    A class for managing sub-cycling times

SourceFiles
    subCycleTime.C

\*---------------------------------------------------------------------------*/

#ifndef subCycleTime_H
#define subCycleTime_H

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class subCycleTime Declaration
\*---------------------------------------------------------------------------*/

class subCycleTime
{
    // Private data

        Time& time_;

        label nSubCycles_;
        label subCycleIndex_;

public:

    // Constructors

        //- Construct from original time and number of sub-cycles
        subCycleTime(Time&, const label nSubCycles);


    //- Destructor
    ~subCycleTime();


    // Member functions

        //- Return true if the number of sub-cycles has been reached
        bool end() const;

        //- End the sub-cycling and reset the time-state
        void endSubCycle();

        //- Return the total number of sub-cycles
        label nSubCycles() const
        {
            return nSubCycles_;
        }

        //- Return the current sub-cycle index
        label index() const
        {
            return subCycleIndex_;
        }


    // Member operators

        //- Prefix increment
        subCycleTime& operator++();

        //- Postfix increment
        subCycleTime& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

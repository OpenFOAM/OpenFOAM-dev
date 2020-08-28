/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::lduScheduleEntry

Description
    Struct to hold the patch index and the initialisation flag for
    the patch schedule

\*---------------------------------------------------------------------------*/

#ifndef lduSchedule_H
#define lduSchedule_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

struct lduScheduleEntry;
Ostream& operator<<(Ostream& os, const lduScheduleEntry& lb);


/*---------------------------------------------------------------------------*\
                         Class lduSchedule Declaration
\*---------------------------------------------------------------------------*/

struct lduScheduleEntry
{
    label patch;
    bool init;

    friend bool operator!=
    (
        const lduScheduleEntry& lb1,
        const lduScheduleEntry& lb2
    )
    {
        return true;
    }

    friend Ostream& operator<<(Ostream& os, const lduScheduleEntry& lb)
    {
        os  << lb.patch << token::SPACE << lb.init;
        return os;
    }
};


typedef List<lduScheduleEntry> lduSchedule;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

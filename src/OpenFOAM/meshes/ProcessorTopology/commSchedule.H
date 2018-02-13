/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::commSchedule

Description
    Determines the order in which a set of processors should communicate
    with one another.

    The communication order should
      - have maximum overlap
      - allow blocking communication without deadlock

    Does a very simple scheduling which assumes same time for all operations.

    After construction:
      - schedule() gives the order in which the input communication should occur
      - procSchedule()[proci] gives per proci

    Does not care whether 'talking' is first send, second receive or maybe
    full swap. This is all responsibility of caller. See ProcessorTopology
    class for use in scheduling processor boundary swaps.

SourceFiles
    commSchedule.C

\*---------------------------------------------------------------------------*/

#ifndef commSchedule_H
#define commSchedule_H

#include "DynamicList.H"
#include "labelPair.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class commSchedule Declaration
\*---------------------------------------------------------------------------*/

class commSchedule
{
    // Private data

        //- Order in which input communication has been scheduled
        labelList schedule_;

        //- Per processor the order in which communication has been scheduled
        labelListList procSchedule_;

    // Private Member Functions

        //- Count the number of outstanding communications for a single
        //  processor
        label outstandingComms(const labelList&, DynamicList<label>&) const;


public:

    ClassName("commSchedule");

    // Constructors

        //- Construct from wanted communication. Wanted communication is between
        //  two processors. Can be a one-way communication or
        //  two-way communication, that is up to the caller. This class just
        //  determines an order for it such that any processor is only talking
        //  to one other at a time. After construction:
        //  - schedule is the order in which comms is done.
        //  - procSchedule[proci] is for proci the order in which comms is done.
        commSchedule(const label nProcs, const List<labelPair>& comms);


    // Member Functions

        //- Order in which comms is scheduled
        const labelList& schedule() const
        {
            return schedule_;
        }

        //- Per processor the order in which communication has been scheduled
        const labelListList& procSchedule() const
        {
            return procSchedule_;
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

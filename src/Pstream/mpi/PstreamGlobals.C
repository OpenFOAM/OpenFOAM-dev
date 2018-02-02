/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

MPI_Comm PstreamGlobals::MPI_COMM_FOAM;

// Outstanding non-blocking operations.
//! \cond fileScope
DynamicList<MPI_Request> PstreamGlobals::outstandingRequests_;
//! \endcond

//// Max outstanding non-blocking operations.
////! \cond fileScope
//int PstreamGlobals::nRequests_ = 0;
////! \endcond

// Free'd non-blocking operations.
//! \cond fileScope
//DynamicList<label> PstreamGlobals::freedRequests_;
//! \endcond

// Max outstanding message tag operations.
//! \cond fileScope
int PstreamGlobals::nTags_ = 0;
//! \endcond

// Free'd message tags
//! \cond fileScope
DynamicList<int> PstreamGlobals::freedTags_;
//! \endcond


// Allocated communicators.
//! \cond fileScope
DynamicList<MPI_Comm> PstreamGlobals::MPICommunicators_;
DynamicList<MPI_Group> PstreamGlobals::MPIGroups_;
//! \endcond

void PstreamGlobals::checkCommunicator
(
    const label comm,
    const label otherProcNo
)
{
    if
    (
        comm < 0
     || comm >= PstreamGlobals::MPICommunicators_.size()
    )
    {
        FatalErrorInFunction
            << "otherProcNo:" << otherProcNo << " : illegal communicator "
            << comm << endl
            << "Communicator should be within range 0.."
            << PstreamGlobals::MPICommunicators_.size()-1 << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

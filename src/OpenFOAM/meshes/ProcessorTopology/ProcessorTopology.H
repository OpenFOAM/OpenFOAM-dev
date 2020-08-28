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
    Foam::ProcessorTopology

Description
    Determines processor-processor connection. After instantiation contains
    on all processors the processor-processor connection table.

    *this[proci] gives the list of neighbouring processors.

    TODO: This does not currently correctly support multiple processor
    patches connecting two processors.

SourceFiles
    ProcessorTopology.C

\*---------------------------------------------------------------------------*/

#ifndef ProcessorTopology_H
#define ProcessorTopology_H

#include "labelList.H"
#include "lduSchedule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ProcessorTopology Declaration
\*---------------------------------------------------------------------------*/

template<class Container, class ProcPatch>
class ProcessorTopology
:
    public labelListList
{
    // Private Data

        //- Local map from neighbour proc to patchi. Different per processor!
        //  -1 or patchi for connection to procID
        labelList procPatchMap_;

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        lduSchedule patchSchedule_;


    // Private Member Functions

        //- Return all neighbouring processors of this processor. Set
        //  procPatchMap_.
        labelList procNeighbours(const label nProcs, const Container&);

public:

    // Constructors

        //- Construct from boundaryMesh
        ProcessorTopology(const Container& patches, const label comm);


    // Member Functions

        //- From neighbour processor to index in boundaryMesh. Local information
        //  (so not same over all processors)
        const labelList& procPatchMap() const
        {
            return procPatchMap_;
        }

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }

        //- Calculate non-blocking (i.e. unscheduled) schedule
        static lduSchedule nonBlockingSchedule(const Container& patches);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ProcessorTopology.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

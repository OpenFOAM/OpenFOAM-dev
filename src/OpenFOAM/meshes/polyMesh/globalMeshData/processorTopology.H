/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

    TODO: This does not currently correctly support multiple processor
    patches connecting two processors.

SourceFiles
    ProcessorTopology.C

\*---------------------------------------------------------------------------*/

#ifndef processorTopology_H
#define processorTopology_H

#include "labelList.H"
#include "lduSchedule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                      Class processorTopology Declaration
\*---------------------------------------------------------------------------*/

class processorTopology
{
    // Private Data

        //- Processor-processor connection table
        labelListList procNbrProcs_;

        //- Local map from neighbour proc to patchi. Different per processor!
        //  -1 or patchi for connection to procID
        labelList procPatchMap_;

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        lduSchedule patchSchedule_;


    // Private Member Functions

        //- Return all neighbouring processors of this processor. Set
        //  procPatchMap_.
        labelList procNeighbours(const label nProcs, const polyBoundaryMesh&);

        //- Calculate non-blocking (i.e. unscheduled) schedule
        lduSchedule nonBlockingSchedule(const polyBoundaryMesh& patches);


public:

    // Constructors

        //- Construct from boundaryMesh
        processorTopology(const polyBoundaryMesh& patches, const label comm);


    // Member Functions

        //- Return the processor-processor connection table
        const labelListList& procNbrProcs() const
        {
            return procNbrProcs_;
        }

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::decompositionMethods::parMetis

Description
    ParMetis redistribution in parallel

    Note: parMetis methods do not support serial operation.

    Parameters
    - Method of decomposition
      - kWay: multilevel k-way
      - geomKway: combined coordinate-based and multi-level k-way
      - adaptiveRepart: balances the work load of a graph

    - Options
      - options[0]: The specified options are used if options[0] = 1

      - options[1]: Speciﬁes the level of information to be returned during
        the execution of the algorithm. Timing information can be obtained by
        setting this to 1. Additional options for this parameter can be obtained
        by looking at parmetis.h. Default: 0.

      - options[2]: Random number seed for the routine

      - options[3]: Speciﬁes whether the sub-domains and processors are coupled
        or un-coupled.  If the number of sub-domains desired (i.e., nparts) and
        the number of processors that are being used is not the same, then these
        must be un-coupled. However, if nparts equals the number of processors,
        these can either be coupled or de-coupled. If sub-domains and processors
        are coupled, then the initial partitioning will be obtained implicitly
        from the graph distribution. However, if sub-domains are un-coupled from
        processors, then the initial partitioning needs to be obtained from the
        initial values assigned to the part array.

    - itr: Parameter which describes the ratio of inter-processor communication
      time compared to data redistribution time.  Should be set between 0.000001
      and 1000000.0.  If set high, a repartitioning with a low edge-cut will be
      computed. If it is set low, a repartitioning that requires little data
      redistribution will be computed.  Good values for this parameter can be
      obtained by dividing inter-processor communication time by data
      redistribution time. Otherwise, a value of 1000.0 is recommended.
      Default: 1000.

SourceFiles
    parMetis.C

\*---------------------------------------------------------------------------*/

#ifndef parMetis_H
#define parMetis_H

#include "decompositionMethod.H"

namespace Foam
{
namespace decompositionMethods
{

/*---------------------------------------------------------------------------*\
                           Class parMetis Declaration
\*---------------------------------------------------------------------------*/

class parMetis
:
    public decompositionMethod
{
    // Private Member Data

        //- Method of decomposition
        word method_;

        //- Options to control the operation of the decomposer
        labelList options_;

        // Processor weights for each constraint
        scalarList processorWeights_;

        //- Parameter which describes the ratio of inter-processor
        //  communication time compared to data redistribution time.
        scalar itr_;


    // Private Member Functions

        label decompose
        (
            const labelList& xadj,
            const labelList& adjncy,
            const pointField& cellCentres,
            const labelList& cellWeights,
            const labelList& faceWeights,
            labelList& finalDecomp
        );


public:

    //- Runtime type information
    TypeName("parMetis");


    // Constructors

        //- Construct given the decomposition dictionary
        parMetis(const dictionary& decompositionDict);

        //- Disallow default bitwise copy construction
        parMetis(const parMetis&) = delete;


    //- Destructor
    virtual ~parMetis()
    {}


    // Member Functions

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        //  Weights get normalised so the minimum value is 1 before truncation
        //  to an integer so the weights should be multiples of the minimum
        //  value. The overall sum of weights might otherwise overflow.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively.
        //  See note on weights above.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const labelList& cellToRegion,
            const pointField& regionPoints,
            const scalarField& regionWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const parMetis&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

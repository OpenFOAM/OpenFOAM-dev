/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::decompositionMethods::ptscotch

Description
    PTScotch domain decomposition.
    For the main details about how to define the strategies, see scotchDecomp.

    Nonetheless, when decomposing in parallel, using <tt>writeGraph=true</tt>
    will write out \c .dgr files for debugging. For example, use these files
    with \c dgpart as follows:

    \verbatim
    mpirun -np 4 dgpart 2 'region0_%r.dgr'
    \endverbatim

    where:
      - %r gets replaced by current processor rank
      - it will decompose into 2 domains

See also
    Foam::decompositionMethods::scotch

SourceFiles
    ptscotch.C

\*---------------------------------------------------------------------------*/

#ifndef ptscotch_H
#define ptscotch_H

#include "decompositionMethod.H"

namespace Foam
{
namespace decompositionMethods
{

/*---------------------------------------------------------------------------*\
                       Class ptscotch Declaration
\*---------------------------------------------------------------------------*/

class ptscotch
:
    public decompositionMethod
{
    dictionary methodDict_;

    // Private Member Functions

        //- Check and print error message
        static void check(const int, const char*);

        //- Decompose. Handles size 0 arrays
        label decompose
        (
            const fileName& meshPath,
            const labelList& adjncy,
            const labelList& xadj,
            const scalarField& cellWeights,
            labelList& decomp
        ) const;

        //- Low level decompose
        label decompose
        (
            const fileName& meshPath,
            const label adjncySize,
            const label adjncy[],
            const label xadjSize,
            const label xadj[],
            const scalarField& cellWeights,
            labelList& decomp
        ) const;


public:

    //- Runtime type information
    TypeName("ptscotch");


    // Constructors

        //- Construct given the decomposition dictionary
        ptscotch
        (
            const dictionary& decompositionDict,
            const dictionary& methodDict
        );

        //- Disallow default bitwise copy construction
        ptscotch(const ptscotch&) = delete;


    //- Destructor
    virtual ~ptscotch()
    {}


    // Member Functions

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed). See note on weights in scotch.H
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively. See note on weights in scotch.H
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const labelList& agglom,
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
        //  See note on weights in scotch.H
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ptscotch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

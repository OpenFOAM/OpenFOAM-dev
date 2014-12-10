/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::scotchDecomp

Description
    Scotch domain decomposition.
    When run in parallel will collect the whole graph on to the master,
    decompose and send back. Run ptscotchDecomp for proper distributed
    decomposition.

SourceFiles
    scotchDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef scotchDecomp_H
#define scotchDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scotchDecomp Declaration
\*---------------------------------------------------------------------------*/

class scotchDecomp
:
    public decompositionMethod
{
    // Private Member Functions

        //- Check and print error message
        static void check(const int, const char*);

        label decompose
        (
            const fileName& meshPath,
            const List<int>& adjncy,
            const List<int>& xadj,
            const scalarField& cWeights,
            List<int>& finalDecomp
        );

        //- Decompose non-parallel
        label decomposeOneProc
        (
            const fileName& meshPath,
            const List<int>& adjncy,
            const List<int>& xadj,
            const scalarField& cWeights,
            List<int>& finalDecomp
        );

        //- Disallow default bitwise copy construct and assignment
        void operator=(const scotchDecomp&);
        scotchDecomp(const scotchDecomp&);


public:

    //- Runtime type information
    TypeName("scotch");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        scotchDecomp(const dictionary& decompositionDict);


    //- Destructor
    virtual ~scotchDecomp()
    {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // Knows about coupled boundaries
            return true;
        }

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        //  Weights get normalised with minimum weight and truncated to
        //  convert into integer so e.g. 3.5 is seen as 3. The overall sum
        //  of weights might otherwise overflow.
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
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

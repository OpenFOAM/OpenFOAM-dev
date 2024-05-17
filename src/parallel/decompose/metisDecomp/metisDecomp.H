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
    Foam::metisDecomp

Description
    Metis domain decomposition

SourceFiles
    metisDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef metisDecomp_H
#define metisDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class metisDecomp Declaration
\*---------------------------------------------------------------------------*/

class metisDecomp
:
    public decompositionMethod
{
    // Private Member Functions

        //- Call Metis with options from dictionary.
        label decompose
        (
            const List<label>& adjncy,
            const List<label>& xadj,
            const scalarField& cellWeights,
            List<label>& finalDecomp
        );


public:

    //- Runtime type information
    TypeName("metis");


    // Constructors

        //- Construct given the decomposition dictionary
        metisDecomp(const dictionary&);

        //- Disallow default bitwise copy construction
        metisDecomp(const metisDecomp&) = delete;


    //- Destructor
    virtual ~metisDecomp()
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
        void operator=(const metisDecomp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

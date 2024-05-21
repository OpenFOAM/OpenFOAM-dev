/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::zoltanDecomp

Description
    Zoltan redistribution in parallel

    Note: Zoltan methods do not support serial operation.

    Parameters
    - lb_method : The load-balancing algorithm
      - block : block partitioning
      - random : random partitioning
      - rcb : recursive coordinate bisection
      - rib : ecursive inertial bisection
      - hsfc : Hilbert space-filling curve partitioning
      - reftree : refinement tree based partitioning
      - graph : choose from collection of methods for graphs
      - hypergraph : choose from a collection of methods for hypergraphs

    - lb_approach The desired load balancing approach. Only lb_method =
      hypergraph or graph uses the lb_approach parameter. Valid values are

      - partition : Partition without reference to the current distribution,
        recommended for static load balancing.

      - repartition : Partition starting from the current data distribution
        to keep data migration low, recommended for dynamic load balancing.

      - refine : Quickly improve the current data distribution

      Default values
      - debug_level     0
      - imbalance_tol   1.05
      - lb_method       graph
      - lb_approach     repartition

Usage
    To select the Zoltan graph repartition method add the following entries to
    decomposeParDict:

        distributor     zoltan;
        libs            ("libzoltanRenumber.so");

    The Zoltan lb_method and lb_approach can be changed by adding the
    corresponding entries to the optional zoltanCeoffs sub-dictionary, e.g.:

    zoltanCoeffs
    {
        lb_method       hypergraph;
        lb_approach     partition;
    }

SourceFiles
    zoltanDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef zoltanDecomp_H
#define zoltanDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class zoltanDecomp Declaration
\*---------------------------------------------------------------------------*/

class zoltanDecomp
:
    public decompositionMethod
{
    // Private Member Functions

        //- Low level decompose
        label decompose
        (
            const CompactListList<label>& adjacency,
            const pointField& points,
            const scalarField& cellWeights,
            List<label>& finalDecomp
        ) const;


public:

    //- Runtime type information
    TypeName("zoltan");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        zoltanDecomp(const dictionary& decompositionDict);

        //- Disallow default bitwise copy construction
        zoltanDecomp(const zoltanDecomp&) = delete;


    //- Destructor
    virtual ~zoltanDecomp()
    {}


    // Member Functions

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed). See note on weights in scotchDecomp.H
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively. See note on weights in scotchDecomp.H
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
        //  See note on weights in scotchDecomp.H
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const zoltanDecomp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

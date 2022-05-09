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
    Foam::nonConformalProcessorCyclicFvPatch

Description
    Non-conformal processor cyclic FV patch. As nonConformalCyclicFvPatch,
    but the neighbouring patch is on a different processor.

See also
    Foam::nonConformalProcessorCyclicPolyPatch
    Foam::nonConformalCyclicFvPatch

SourceFiles
    nonConformalProcessorCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicFvPatch_H
#define nonConformalProcessorCyclicFvPatch_H

#include "processorCyclicFvPatch.H"
#include "nonConformalProcessorCyclicPolyPatch.H"
#include "nonConformalCoupledFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class nonConformalProcessorCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalProcessorCyclicFvPatch
:
    public processorCyclicFvPatch,
    public nonConformalCoupledFvPatch
{
    // Private Data

        //- Reference to the polyPatch
        const nonConformalProcessorCyclicPolyPatch&
            nonConformalProcessorCyclicPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField& w) const;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        nonConformalProcessorCyclicFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~nonConformalProcessorCyclicFvPatch();


    // Member Functions

        // Access

            //- Poly patch
            const nonConformalProcessorCyclicPolyPatch&
                nonConformalProcessorCyclicPatch() const;

            //- ...
            using processorCyclicFvPatch::owner;

            //- ...
            using processorCyclicFvPatch::neighbour;

            //- ...
            using processorCyclicFvPatch::transform;

            //- Is this patch coupled? Returns true. For NCC patches the poly
            //  mesh is considered non-coupled whilst the finite volume mesh is
            //  considered coupled.
            virtual bool coupled() const;

            //- Return the start label of this patch in the polyMesh face list.
            //  Raises an error unless the patch is empty.
            virtual label start() const;

            //- Return the size
            virtual label size() const;

            //- Return the face-cells
            virtual const labelUList& faceCells() const;


        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

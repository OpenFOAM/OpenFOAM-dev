/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::nonConformalCyclicLagrangianPatch

Description
    Non-conformal-cyclic Lagrangian patch

SourceFiles
    nonConformalCyclicLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicLagrangianPatch_H
#define nonConformalCyclicLagrangianPatch_H

#include "LagrangianPatch.H"
#include "nonConformalCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nonConformalCyclicLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCyclicLagrangianPatch
:
    public LagrangianPatch
{
private:

    // Private Data

        //- Reference to the nonConformalCyclic poly patch
        const nonConformalCyclicPolyPatch& nonConformalCyclicPatch_;

        //- Between crossFaces/evaluate and partition the elements have moved
        //  to the neighbour patch, but the mesh has not been partitioned. So,
        //  the sub-mesh that relates to a cyclic patch's elements is that of
        //  the neighbour patch. This flag identifies when this is the case.
        mutable bool isNbrPatchMesh_;

        //- The number of boundary hits that occurred during locate executions
        mutable label nPositionalErrors_;

        //- The maximum positional error (squared) following a transfer
        mutable scalar maxPositionalErrorSqr_;

        //- The location of the maximum positional error following a transfer
        mutable point maxPositionalErrorReceivePosition_;


public:

    // Friend classes

        //- The non-conformal processor cyclic needs to add into the positional
        //  error variables of the non-processor cyclic patch
        friend class nonConformalProcessorCyclicLagrangianPatch;


    //- Runtime type information
    TypeName(nonConformalCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        nonConformalCyclicLagrangianPatch
        (
            const polyPatch&,
            const LagrangianBoundaryMesh&
        );


    //- Destructor
    virtual ~nonConformalCyclicLagrangianPatch();


    // Member Functions

        //- Return the nonConformalCyclic patch
        const nonConformalCyclicPolyPatch& nonConformalCyclicPatch() const
        {
            return nonConformalCyclicPatch_;
        }

        //- Return a reference to the nonConformalCyclic patch's transform
        const transformer& transform() const
        {
            return nonConformalCyclicPatch_.transform();
        }

        //- Return the sub-mesh associated with this patch
        virtual const LagrangianSubMesh& mesh() const;

        //- Evaluate changes in elements that have tracked to this patch
        virtual void evaluate
        (
            PstreamBuffers&,
            LagrangianMesh&,
            const LagrangianScalarInternalDynamicField& fraction
        ) const;

        //- Update following partitioning of the mesh
        virtual void partition() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

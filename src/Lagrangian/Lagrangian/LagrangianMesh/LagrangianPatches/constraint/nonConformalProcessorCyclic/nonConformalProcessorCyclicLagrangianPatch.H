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
    Foam::nonConformalProcessorCyclicLagrangianPatch

Description
    Non-conformal-processor-cyclic Lagrangian patch

SourceFiles
    nonConformalProcessorCyclicLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicLagrangianPatch_H
#define nonConformalProcessorCyclicLagrangianPatch_H

#include "processorCyclicLagrangianPatch.H"
#include "nonConformalCyclicLagrangianPatch.H"
#include "nonConformalProcessorCyclicPolyPatch.H"
#include "LagrangianBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nonConformalProcessorCyclicLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalProcessorCyclicLagrangianPatch
:
    public processorCyclicLagrangianPatch
{
private:

    // Private Data

        //- Reference to the nonConformalProcessorCyclic poly patch
        const nonConformalProcessorCyclicPolyPatch&
            nonConformalProcessorCyclicPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        nonConformalProcessorCyclicLagrangianPatch
        (
            const polyPatch&,
            const LagrangianBoundaryMesh&
        );


    //- Destructor
    virtual ~nonConformalProcessorCyclicLagrangianPatch();


    // Member Functions

        //- Return the nonConformalProcessorCyclic patch
        const nonConformalProcessorCyclicPolyPatch&
        nonConformalProcessorCyclicPatch() const
        {
            return nonConformalProcessorCyclicPatch_;
        }

        //- Return a reference to the referred patch
        const nonConformalCyclicLagrangianPatch& referPatch() const
        {
            return
                refCast<const nonConformalCyclicLagrangianPatch>
                (
                    boundaryMesh()[referPatchIndex()]
                );
        }

        //- Initialise evaluation of changes in elements that have tracked to
        //  this patch
        virtual void initEvaluate
        (
            PstreamBuffers&,
            LagrangianMesh&,
            const LagrangianScalarInternalDynamicField& fraction
        ) const;

        //- Evaluate changes in elements that have tracked to this patch
        virtual void evaluate
        (
            PstreamBuffers&,
            LagrangianMesh&,
            const LagrangianScalarInternalDynamicField& fraction
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

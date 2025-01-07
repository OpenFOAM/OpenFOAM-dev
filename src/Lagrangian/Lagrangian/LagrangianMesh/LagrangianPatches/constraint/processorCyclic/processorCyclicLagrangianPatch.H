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
    Foam::processorCyclicLagrangianPatch

Description
    Processor-cyclic Lagrangian patch. This is used for the patches that
    interface between processors across cyclic patches in the complete mesh.

SourceFiles
    processorCyclicLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicLagrangianPatch_H
#define processorCyclicLagrangianPatch_H

#include "processorLagrangianPatch.H"
#include "cyclicLagrangianPatch.H"
#include "processorCyclicPolyPatch.H"
#include "LagrangianBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorCyclicLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicLagrangianPatch
:
    public processorLagrangianPatch
{
private:

    // Private Data

        //- Reference to the processorCyclic poly patch
        const processorCyclicPolyPatch& processorCyclicPatch_;


public:

    //- Runtime type information
    TypeName(processorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        processorCyclicLagrangianPatch
        (
            const polyPatch&,
            const LagrangianBoundaryMesh&
        );


    //- Destructor
    virtual ~processorCyclicLagrangianPatch();


    // Member Functions

        //- Return the processorCyclic patch
        const processorCyclicPolyPatch& processorCyclicPatch() const
        {
            return processorCyclicPatch_;
        }

        //- Return the referred patch index
        label referPatchIndex() const
        {
            return processorCyclicPatch_.referPatchIndex();
        }

        //- Return a reference to the referred patch
        const cyclicLagrangianPatch& referPatch() const
        {
            return
                refCast<const cyclicLagrangianPatch>
                (
                    boundaryMesh()[referPatchIndex()]
                );
        }

        //- Return a reference to the referred patch's transform
        const transformer& transform() const
        {
            return processorCyclicPatch_.transform();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

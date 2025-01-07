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
    Foam::processorLagrangianPatch

Description
    Processor Lagrangian patch. This is used for the patches that interface
    between processors across faces that are internal to the complete mesh.

    with a co-located patch on another process.

SourceFiles
    processorLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorLagrangianPatch_H
#define processorLagrangianPatch_H

#include "LagrangianPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class processorLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class processorLagrangianPatch
:
    public LagrangianPatch
{
private:

    // Private Data

        //- Reference to the processor poly patch
        const processorPolyPatch& processorPatch_;


protected:

    // Protected Data

        //- The sub-mesh used to append received data
        mutable autoPtr<LagrangianSubMesh> receiveMeshPtr_;


public:

    //- Runtime type information
    TypeName(processorPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        processorLagrangianPatch
        (
            const polyPatch&,
            const LagrangianBoundaryMesh&
        );


    //- Destructor
    virtual ~processorLagrangianPatch();


    // Member Functions

        //- Return the processor patch
        const processorPolyPatch& processorPatch() const
        {
            return processorPatch_;
        }

        //- Return the sub-mesh associated with this patch
        virtual const LagrangianSubMesh& mesh() const;

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

        //- Update following partitioning of the mesh
        virtual void partition() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

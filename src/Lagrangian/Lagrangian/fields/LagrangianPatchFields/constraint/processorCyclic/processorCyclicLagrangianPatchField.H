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
    Foam::processorCyclicLagrangianPatchField

Description
    A processorCyclic boundary condition for Lagrangian. Properties are
    communicated to and from the neighbour patch and are transformed by the
    transformation of the associated cyclic patch.

Usage
    Example specification:
    <LagrangianPatchName>
    {
        type            processorCyclic;
    }
    \endverbatim

SourceFiles
    processorCyclicLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicLagrangianPatchField_H
#define processorCyclicLagrangianPatchField_H

#include "processorLagrangianPatchField.H"
#include "processorCyclicLagrangianPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class processorCyclicLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorCyclicLagrangianPatchField
:
    public processorLagrangianPatchField<Type>
{
private:

    // Private Data

        //- Reference to the processorCyclic patch
        const processorCyclicLagrangianPatch& processorCyclicPatch_;


public:

    //- Runtime type information
    TypeName("processorCyclic");


    // Constructors

        //- Construct from patch and internal field
        processorCyclicLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        processorCyclicLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        processorCyclicLagrangianPatchField
        (
            const processorCyclicLagrangianPatchField<Type>&
        );

        //- Construct as copy setting the internal field reference
        processorCyclicLagrangianPatchField
        (
            const processorCyclicLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct as copy setting the patch and the internal field reference
        processorCyclicLagrangianPatchField
        (
            const processorCyclicLagrangianPatchField<Type>&,
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new processorCyclicLagrangianPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianPatchField<Type>> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new processorCyclicLagrangianPatchField<Type>
                (
                    *this,
                    iIo
                )
            );
        }

        //- Construct and return a clone onto a new patch
        //  and setting the internal field reference
        virtual autoPtr<LagrangianPatchField<Type>> clone
        (
            const LagrangianPatch& p,
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new processorCyclicLagrangianPatchField<Type>
                (
                    *this,
                    p,
                    iIo
                )
            );
        }


    // Member Functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            PstreamBuffers&,
            const LagrangianScalarInternalDynamicField& fraction
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorCyclicLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::processorLagrangianPatchField

Description
    A processor boundary condition for Lagrangian. Properties are communicated
    to and from the neighbour patch.

Usage
    Example specification:
    <LagrangianPatchName>
    {
        type            processor;
    }
    \endverbatim

SourceFiles
    processorLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef processorLagrangianPatchField_H
#define processorLagrangianPatchField_H

#include "LagrangianPatchField.H"
#include "processorLagrangianPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class processorLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorLagrangianPatchField
:
    public LagrangianPatchField<Type>
{
private:

    // Private Data

        //- Reference to the processor patch
        const processorLagrangianPatch& processorPatch_;


    // Private Member Functions

        //- Evaluate, given the internal field
        template<template<class> class GeoField>
        void evaluate(PstreamBuffers&, const GeoField<Type>& internalField);


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from patch and internal field
        processorLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        processorLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        processorLagrangianPatchField
        (
            const processorLagrangianPatchField<Type>&
        );

        //- Construct as copy setting the internal field reference
        processorLagrangianPatchField
        (
            const processorLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct as copy setting the patch and the internal field reference
        processorLagrangianPatchField
        (
            const processorLagrangianPatchField<Type>&,
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new processorLagrangianPatchField<Type>
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
                new processorLagrangianPatchField<Type>
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
                new processorLagrangianPatchField<Type>
                (
                    *this,
                    p,
                    iIo
                )
            );
        }


    // Member Functions

        //- Initialise evaluation of the patch field
        virtual void initEvaluate
        (
            PstreamBuffers&,
            const LagrangianScalarInternalDynamicField& fraction
        );

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
    #include "processorLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::nonConformalProcessorCyclicLagrangianPatchField

Description
    A non-conformal cyclic boundary condition for Lagrangian. Properties are
    transformed by the transformation of the associated non-conformal cyclic
    patch.

Usage
    Example specification:
    <LagrangianPatchName>
    {
        type            nonConformalProcessorCyclic;
    }
    \endverbatim

SourceFiles
    nonConformalProcessorCyclicLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicLagrangianPatchField_H
#define nonConformalProcessorCyclicLagrangianPatchField_H

#include "processorCyclicLagrangianPatchField.H"
#include "nonConformalProcessorCyclicLagrangianPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class nonConformalProcessorCyclicLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalProcessorCyclicLagrangianPatchField
:
    public processorCyclicLagrangianPatchField<Type>
{
private:

    // Private Data

        //- Reference to the nonConformalProcessorCyclic patch
        const nonConformalProcessorCyclicLagrangianPatch&
            nonConformalProcessorCyclicPatch_;


public:

    //- Runtime type information
    TypeName("nonConformalProcessorCyclic");


    // Constructors

        //- Construct from patch and internal field
        nonConformalProcessorCyclicLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalProcessorCyclicLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        nonConformalProcessorCyclicLagrangianPatchField
        (
            const nonConformalProcessorCyclicLagrangianPatchField<Type>&
        );

        //- Construct as copy setting the internal field reference
        nonConformalProcessorCyclicLagrangianPatchField
        (
            const nonConformalProcessorCyclicLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct as copy setting the patch and the internal field reference
        nonConformalProcessorCyclicLagrangianPatchField
        (
            const nonConformalProcessorCyclicLagrangianPatchField<Type>&,
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new nonConformalProcessorCyclicLagrangianPatchField<Type>
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
                new nonConformalProcessorCyclicLagrangianPatchField<Type>
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
                new nonConformalProcessorCyclicLagrangianPatchField<Type>
                (
                    *this,
                    p,
                    iIo
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalProcessorCyclicLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

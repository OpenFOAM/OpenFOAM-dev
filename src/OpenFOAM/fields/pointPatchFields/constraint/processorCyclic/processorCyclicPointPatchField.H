/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::processorCyclicPointPatchField

Description
    Foam::processorCyclicPointPatchField

SourceFiles
    processorCyclicPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicPointPatchField_H
#define processorCyclicPointPatchField_H

#include "coupledPointPatchField.H"
#include "processorCyclicPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class processorCyclicPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorCyclicPointPatchField
:
    public coupledPointPatchField<Type>
{
    // Private Data

        //- Local reference to processor patch
        const processorCyclicPointPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(processorCyclicPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        processorCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        processorCyclicPointPatchField
        (
            const processorCyclicPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        processorCyclicPointPatchField
        (
            const processorCyclicPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new processorCyclicPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }

        //- Construct and return a clone onto a new patch
        //  setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const pointPatch& patch,
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new processorCyclicPointPatchField<Type>
                (
                    patch,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~processorCyclicPointPatchField();


    // Member Functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const
            {
                if (Pstream::parRun())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}


    // Member Operators

        //- Inherit assignment
        using coupledPointPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorCyclicPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

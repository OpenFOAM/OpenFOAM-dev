/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::processorCyclicFvsPatchField

Description
    Foam::processorCyclicFvsPatchField

SourceFiles
    processorCyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicFvsPatchField_H
#define processorCyclicFvsPatchField_H

#include "processorCyclicFvPatch.H"
#include "processorFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class processorCyclicFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorCyclicFvsPatchField
:
    public processorFvsPatchField<Type>
{
    // Private Data

        //- Local reference cast into the processor patch
        const processorCyclicFvPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(processorCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given processorCyclicFvsPatchField onto a
        //  new patch
        processorCyclicFvsPatchField
        (
            const processorCyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        processorCyclicFvsPatchField
        (
            const processorCyclicFvsPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        processorCyclicFvsPatchField
        (
            const processorCyclicFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new processorCyclicFvsPatchField<Type>(*this, iF)
            );
        }

        //- Construct and return a clone onto a new patch
        //  setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const fvPatch& patch,
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new processorCyclicFvsPatchField<Type>(patch, iF, *this)
            );
        }


    //- Destructor
    virtual ~processorCyclicFvsPatchField();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorCyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

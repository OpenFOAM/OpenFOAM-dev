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
    Foam::nonConformalProcessorCyclicFvsPatchField

Description
    Foam::nonConformalProcessorCyclicFvsPatchField

SourceFiles
    nonConformalProcessorCyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicFvsPatchField_H
#define nonConformalProcessorCyclicFvsPatchField_H

#include "nonConformalProcessorCyclicFvPatch.H"
#include "processorCyclicFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class nonConformalProcessorCyclicFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalProcessorCyclicFvsPatchField
:
    public processorCyclicFvsPatchField<Type>
{
    // Private Data

        //- Local reference cast into the processor patch
        const nonConformalProcessorCyclicFvPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalProcessorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch and internal field and patch field
        nonConformalProcessorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalProcessorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalProcessorCyclicFvsPatchField
        //  onto a new patch
        nonConformalProcessorCyclicFvsPatchField
        (
            const nonConformalProcessorCyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalProcessorCyclicFvsPatchField
        (
            const nonConformalProcessorCyclicFvsPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalProcessorCyclicFvsPatchField
        (
            const nonConformalProcessorCyclicFvsPatchField<Type>&,
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
                new nonConformalProcessorCyclicFvsPatchField<Type>(*this, iF)
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
                new nonConformalProcessorCyclicFvsPatchField<Type>
                (
                    patch,
                    iF,
                    *this
                )
            );
        }


    //- Destructor
    virtual ~nonConformalProcessorCyclicFvsPatchField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalProcessorCyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

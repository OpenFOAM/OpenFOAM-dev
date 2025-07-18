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
    Foam::codedFixedValueFvPatchField

Description
    Constructs on-the-fly a new boundary condition (derived from
    fixedValueFvPatchField) which is then used to evaluate.

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            codedFixedValue;
        value           uniform 0;
        name    rampedFixedValue;   // name of generated BC

        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};

        // codeInclude
        //#{
        //    #include ".H"
        //#};

        // codeOptions
        //#{
        //    -I$(LIB_SRC)/finiteVolume/lnInclude
        //#};
    }
    \endverbatim

    A special form is if the 'code' section is not supplied. In this case
    the code is read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry:

    \verbatim
    <patchName>
    {
        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};
    }
    \endverbatim

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

SourceFiles
    codedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef codedFixedValueFvPatchField_H
#define codedFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class codedFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public codedBase
{
    // Private Data

        //- Keywords associated with source code
        static const wordList codeKeys;

        //- Name of the dictionary variables in the source code
        static const wordList codeDictVars;

        mutable autoPtr<fvPatchField<Type>> redirectPatchFieldPtr_;


    // Private Member Functions

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;


public:

    //- Runtime type information
    TypeName("codedFixedValue");


    // Constructors

        //- Construct from patch, internal field and dictionary
        codedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given codedFixedValueFvPatchField
        //  onto a new patch
        codedFixedValueFvPatchField
        (
            const codedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );


        //- Disallow copy without setting internal field reference
        codedFixedValueFvPatchField
        (
            const codedFixedValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        codedFixedValueFvPatchField
        (
            const codedFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new codedFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Get reference to the underlying patch
        const fvPatchField<Type>& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field, sets Updated to false
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "codedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

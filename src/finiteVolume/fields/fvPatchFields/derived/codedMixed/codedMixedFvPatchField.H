/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::codedMixedFvPatchField

Description
    Mixed boundary condition evaluated from the code provided

    which is compiled, linked and executed automatically at run-time.

    The optional \c name entry specifies the name of the code generated from the
    template, if not provided the it constructed as
    \<field name\>_\<patch name\> for single region cases or
    \<region name\>_\<field name\>_\<patch name\> for multi-region cases.

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            codedMixed;

        refValue        uniform (0 0 0);
        refGradient     uniform (0 0 0);
        valueFraction   uniform 1;

        name            rampedFixedValue;   // Optional name of generated code

        code
        #{
            this->refValue() =
                vector(1, 0, 0)*min(10, 0.1*this->time().value());
            this->refGrad() = Zero;
            this->valueFraction() = 1.0;
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

See also
    Foam::codedBase
    Foam::dynamicCode

SourceFiles
    codedMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef codedMixedFvPatchField_H
#define codedMixedFvPatchField_H

#include "mixedFvPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class codedMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedMixedFvPatchField
:
    public mixedFvPatchField<Type>,
    public codedBase
{
    // Private Data

        //- Keywords associated with source code
        static const wordList codeKeys;

        //- Name of the dictionary variables in the source code
        static const wordList codeDictVars;

        // Name of the code options template file
        static const word codeOptions;

        // Names of the source files to compile
        static const wordList compileFiles;

        // Names of the source files to copy (header files)
        static const wordList copyFiles;

        mutable autoPtr<mixedFvPatchField<Type>> redirectPatchFieldPtr_;


public:

    //- Runtime type information
    TypeName("codedMixed");


    // Constructors

        //- Construct from patch, internal field and dictionary
        codedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given codedMixedFvPatchField
        //  onto a new patch
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const fieldMapper&
        );


        //- Disallow copy without setting internal field reference
        codedMixedFvPatchField(const codedMixedFvPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&,
            const DimensionedField<Type, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new codedMixedFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        using fvPatchField<Type>::patch;
        using fvPatchField<Type>::internalField;


        //- Get reference to the underlying patchField
        const mixedFvPatchField<Type>& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        //  This is only needed to set the updated() flag of the name
        //  to false.
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using mixedFvPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "codedMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

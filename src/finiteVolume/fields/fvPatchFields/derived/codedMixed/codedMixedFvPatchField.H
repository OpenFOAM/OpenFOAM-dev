/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Constructs on-the-fly a new boundary condition (derived from
    mixedFvPatchField) which is then used to evaluate.

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            codedMixed;

        refValue        uniform (0 0 0);
        refGradient     uniform (0 0 0);
        valueFraction   uniform 1;

        name    rampedMixed;   // name of generated BC

        code
        #{
            this->refValue() =
                vector(1, 0, 0)
               *min(10, 0.1*this->db().time().value());
            this->refGrad() = Zero;
            this->valueFraction() = 1.0;
        #};

        // codeInclude
        //#{
        //    #include "fvCFD.H"
        //#};

        // codeOptions
        //#{
        //    -I$(LIB_SRC)/finiteVolume/lnInclude
        //#};
    }
    \endverbatim

    A special form is if the 'code' section is not supplied. In this case
    the code gets read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry

    \verbatim
    <patchName>
    {
        code
        #{
            this->refValue() = min(10, 0.1*this->db().time().value());
            this->refGrad() = Zero;
            this->valueFraction() = 1.0;
        #};
    }
    \endverbatim

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

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

// Forward declaration of classes
class dynamicCode;
class dynamicCodeContext;
class IOdictionary;

/*---------------------------------------------------------------------------*\
                   Class codedMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedMixedFvPatchField
:
    public mixedFvPatchField<Type>,
    public codedBase
{
    // Private data

        //- Dictionary contents for the boundary condition
        mutable dictionary dict_;

        const word name_;

        mutable autoPtr<mixedFvPatchField<Type>> redirectPatchFieldPtr_;


    // Private Member Functions

        const IOdictionary& dict() const;

        //- Set the rewrite vars controlling the Type
        static void setFieldTemplates(dynamicCode& dynCode);

        //- Get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        // Return a description (type + name) for the output
        virtual string description() const;

        // Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;


public:

    // Static data members

        //- Name of the C code template to be used
        static const word codeTemplateC;

        //- Name of the H code template to be used
        static const word codeTemplateH;


    //- Runtime type information
    TypeName("codedMixed");


    // Constructors

        //- Construct from patch and internal field
        codedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given codedMixedFvPatchField
        //  onto a new patch
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new codedMixedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&,
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
                new codedMixedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

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

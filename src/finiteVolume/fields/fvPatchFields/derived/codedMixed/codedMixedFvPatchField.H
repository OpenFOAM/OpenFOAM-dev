/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
        //    #include ".H"
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

        mutable autoPtr<mixedFvPatchField<Type>> redirectPatchFieldPtr_;


    // Private Member Functions

        //- Get the keywords associated with source code
        virtual wordList codeKeys() const;

        //- Get the name of the dictionary variables in the source code
        virtual wordList codeDictVars() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const;


public:

    //- Runtime type information
    TypeName("codedMixed");


    // Constructors

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
            const fieldMapper&
        );


        //- Disallow copy without setting internal field reference
        codedMixedFvPatchField(const codedMixedFvPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
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


    // Member Functions

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

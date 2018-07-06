/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::codedFixedValuePointPatchField

Description
    Constructs on-the-fly a new boundary condition (derived from
    fixedValuePointPatchField) which is then used to evaluate.

    Example:
    \verbatim
    movingWall
    {
        type            codedFixedValue;
        value           uniform 0;
        name    rampedFixedValue;   // name of generated bc

        code
        #{
            operator==
            (
                vector(0,0,1)
                *min(10, 0.1*this->db().time().value())
            );
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

    A special form is if the \c code section is not supplied. In this case
    the code gets read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry

    \verbatim
    rampedFixedValue
    {
        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};
    }
    \endverbatim

See also
    codedFixedValueFvPatchField

SourceFiles
    codedFixedValuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef codedFixedValuePointPatchField_H
#define codedFixedValuePointPatchField_H

#include "fixedValuePointPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicCode;
class dynamicCodeContext;
class IOdictionary;

/*---------------------------------------------------------------------------*\
               Class codedFixedValuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedFixedValuePointPatchField
:
    public fixedValuePointPatchField<Type>,
    public codedBase
{
    // Private data

        //- Dictionary contents for the boundary condition
        mutable dictionary dict_;

        const word name_;

        mutable autoPtr<pointPatchField<Type>> redirectPatchFieldPtr_;


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
    TypeName("codedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        codedFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given codedFixedValuePointPatchField
        //  onto a new patch
        codedFixedValuePointPatchField
        (
            const codedFixedValuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        codedFixedValuePointPatchField
        (
            const codedFixedValuePointPatchField<Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new codedFixedValuePointPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        codedFixedValuePointPatchField
        (
            const codedFixedValuePointPatchField<Type>&,
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
                new codedFixedValuePointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get reference to the underlying patch
        const pointPatchField<Type>& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field, sets Updated to false
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
    #include "codedFixedValuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

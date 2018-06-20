/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::outletMappedUniformInletFvPatchField

Description
    This boundary conditon averages the field over the "outlet" patch specified
    by name "outletPatch" and applies this as the uniform value of the field
    over this patch.

Usage
    \table
        Property        | Description             | Required    | Default value
        outletPatch     | Name of outlet patch    | yes         |
        phi             | Flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            outletMappedUniformInlet;
        outletPatchName aPatch;
        phi             phi;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    outletMappedUniformInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef outletMappedUniformInletFvPatchField_H
#define outletMappedUniformInletFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class outletMappedUniformInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class outletMappedUniformInletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Name of the outlet patch to be mapped
        word outletPatchName_;

       //- Name of the flux transporting the field
        word phiName_;


public:

    //- Runtime type information
    TypeName("outletMappedUniformInlet");


    // Constructors

        //- Construct from patch and internal field
        outletMappedUniformInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletMappedUniformInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletMappedUniformInletFvPatchField
        //  onto a new patch
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new outletMappedUniformInletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&,
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
                new outletMappedUniformInletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Name of the outlet patch to be mapped
            const word& outletPatchName() const
            {
                return outletPatchName_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "outletMappedUniformInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

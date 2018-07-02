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
    Foam::slipFvPatchField

Description
    This boundary condition provides a slip constraint.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            slip;
    }
    \endverbatim

SourceFiles
    slipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef slipFvPatchField_H
#define slipFvPatchField_H

#include "basicSymmetryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class slipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class slipFvPatchField
:
    public basicSymmetryFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("slip");


    // Constructors

        //- Construct from patch and internal field
        slipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        slipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given slipFvPatchField onto a new patch
        slipFvPatchField
        (
            const slipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        slipFvPatchField
        (
            const slipFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new slipFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        slipFvPatchField
        (
            const slipFvPatchField<Type>&,
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
                new slipFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


    // Member operators

        // virtual void operator=(const UList<Type>&) {}

        // virtual void operator=(const fvPatchField<Type>&) {}
        // virtual void operator+=(const fvPatchField<Type>&) {}
        // virtual void operator-=(const fvPatchField<Type>&) {}
        // virtual void operator*=(const fvPatchField<scalar>&) {}
        // virtual void operator/=(const fvPatchField<scalar>&) {}

        // virtual void operator+=(const Field<Type>&) {}
        // virtual void operator-=(const Field<Type>&) {}

        // virtual void operator*=(const Field<scalar>&) {}
        // virtual void operator/=(const Field<scalar>&) {}

        // virtual void operator=(const Type&) {}
        // virtual void operator+=(const Type&) {}
        // virtual void operator-=(const Type&) {}
        // virtual void operator*=(const scalar) {}
        // virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "slipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

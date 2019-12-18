/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::mixedFvPatchField

Description
    This boundary condition provides a base class for 'mixed' type boundary
    conditions, i.e. conditions that mix fixed value and patch-normal gradient
    conditions.

    The respective contributions from each is determined by a weight field:

        \f[
            x_p = w x_p + (1-w) \left(x_c + \frac{\nabla_\perp x}{\Delta}\right)
        \f]

    where
    \vartable
        x_p   | patch values
        x_c   | patch internal cell values
        w     | weight field
        \Delta| inverse distance from face centre to internal cell centre
        w     | weighting (0-1)
    \endvartable


Usage
    \table
        Property     | Description             | Required    | Default value
        valueFraction | weight field           | yes         |
        refValue     | fixed value             | yes         |
        refGrad      | patch normal gradient   | yes         |
    \endtable

Note
    This condition is not usually applied directly; instead, use a derived
    mixed condition such as \c inletOutlet

See also
    Foam::inletOutletFvPatchField

SourceFiles
    mixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedFvPatchField_H
#define mixedFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixedFvPatchField
:
    public fvPatchField<Type>
{
    // Private Data

        //- Value field
        Field<Type> refValue_;

        //- Normal gradient field
        Field<Type> refGrad_;

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;


public:

    //- Runtime type information
    TypeName("mixed");


    // Constructors

        //- Construct from patch and internal field
        mixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given mixedFvPatchField onto a new patch
        mixedFvPatchField
        (
            const mixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Copy constructor
        mixedFvPatchField
        (
            const mixedFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mixedFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        mixedFvPatchField
        (
            const mixedFvPatchField<Type>&,
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
                new mixedFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


        // Return defining fields

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }

            virtual Field<Type>& refGrad()
            {
                return refGrad_;
            }

            virtual const Field<Type>& refGrad() const
            {
                return refGrad_;
            }

            virtual scalarField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchField<Type>&, const labelList&);


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

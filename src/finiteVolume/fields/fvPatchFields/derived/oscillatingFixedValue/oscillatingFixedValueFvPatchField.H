/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::oscillatingFixedValueFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an oscillating condition in terms of
    amplitude and frequency.

        \f[
            x_p = (1 + a sin(2 \pi f t))x_{ref} + x_o
        \f]

    where

    \vartable
        x_p     | patch values
        x_{ref} | patch reference values
        x_o     | patch offset values
        a       | amplitude
        f       | frequency [1/s]
        t       | time [s]
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        refValue     | reference value         | yes         |
        offset       | offset value            | no          | 0.0
        amplitude    | oscillation amplitude   | yes         |
        frequency    | oscillation frequency   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            oscillatingFixedValue;
        refValue        uniform 5.0;
        offset          0.0;
        amplitude       constant 0.5;
        frequency       constant 10;
    }
    \endverbatim

Note
    The amplitude and frequency entries are DataEntry types, able to describe
    time varying functions.  The example above gives the usage for supplying
    constant values.

SeeAlso
    Foam::DataEntry

SourceFiles
    oscillatingFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingFixedValueFvPatchField_H
#define oscillatingFixedValueFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class oscillatingFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class oscillatingFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Reference value
        Field<Type> refValue_;

        //- Offset
        Type offset_;

        //- Amplitude
        autoPtr<DataEntry<scalar> > amplitude_;

        //- Frequency
        autoPtr<DataEntry<scalar> > frequency_;

        //- Current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Return current scale
        scalar currentScale() const;


public:

    //- Runtime type information
    TypeName("oscillatingFixedValue");


    // Constructors

        //- Construct from patch and internal field
        oscillatingFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oscillatingFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oscillatingFixedValueFvPatchField
        //  onto a new patch
        oscillatingFixedValueFvPatchField
        (
            const oscillatingFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oscillatingFixedValueFvPatchField
        (
            const oscillatingFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new oscillatingFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        oscillatingFixedValueFvPatchField
        (
            const oscillatingFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new oscillatingFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            const Field<Type>& refValue() const
            {
                return refValue_;
            }

            //- Return reference to the ref value to allow adjustment
            Field<Type>& refValue()
            {
                return refValue_;
            }

            //- Return amplitude
            scalar amplitude() const
            {
                return amplitude_;
            }

            scalar& amplitude()
            {
                return amplitude_;
            }

            //- Return frequency
            scalar frequency() const
            {
                return frequency_;
            }

            scalar& frequency()
            {
                return frequency_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


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
#   include "oscillatingFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

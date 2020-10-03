/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::pressureFvPatchScalarField

Description
    Static pressure boundary condition

    Used as the base-class for static pressure boundary conditions force
    alternative pressure variables, e.g. p_rgh.

Usage
    \table
        Property     | Description             | Required    | Default value
        p            | static pressure         | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressure;
        p               uniform 0;
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchScalarField
    Foam::PrghPressureFvPatchScalarField

SourceFiles
    pressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureFvPatchScalarField_H
#define pressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class pressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class pressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data

        //- Static pressure
        scalarField p_;


public:

    //- Runtime type information
    TypeName("pressure");


    // Constructors

        //- Construct from patch and internal field
        pressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureFvPatchScalarField onto a new patch
        pressureFvPatchScalarField
        (
            const pressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        pressureFvPatchScalarField(const pressureFvPatchScalarField&) = delete;

        //- Copy constructor setting internal field reference
        pressureFvPatchScalarField
        (
            const pressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the static pressure
            const scalarField& p() const
            {
                return p_;
            }

            //- Return reference to the static pressure to allow adjustment
            scalarField& p()
            {
                return p_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

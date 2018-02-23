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
    Foam::freestreamPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a free-stream condition for pressure.

    It is an outlet-inlet condition that uses the velocity orientation to
    continuously blend between zero gradient for normal inlet and fixed value
    for normal outlet flow.

Usage
    \table
        Property        | Description             | Required    | Default value
        U               | velocity field name     | no          | U
        freestreamValue | freestream pressure     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            freestreamPressure;
        freestreamValue uniform 1e5;
    }
    \endverbatim

Note
    This condition is designed to operate with a freestreamVelocity condition

See also
    Foam::mixedFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    freestreamPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef freestreamPressureFvPatchScalarField_H
#define freestreamPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class freestreamPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class freestreamPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("freestreamPressure");


    // Constructors

        //- Construct from patch and internal field
        freestreamPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freestreamPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given freestreamPressureFvPatchScalarField onto
        //  a new patch
        freestreamPressureFvPatchScalarField
        (
            const freestreamPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        freestreamPressureFvPatchScalarField
        (
            const freestreamPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new freestreamPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        freestreamPressureFvPatchScalarField
        (
            const freestreamPressureFvPatchScalarField&,
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
                new freestreamPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

            const scalarField& freestreamValue() const
            {
                return refValue();
            }

            scalarField& freestreamValue()
            {
                return refValue();
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

#endif

// ************************************************************************* //

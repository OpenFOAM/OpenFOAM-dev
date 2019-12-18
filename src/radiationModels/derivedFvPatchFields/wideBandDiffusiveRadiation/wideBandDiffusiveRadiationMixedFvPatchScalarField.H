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
    Foam::wideBandDiffusiveRadiationMixedFvPatchScalarField

Description
    This boundary condition provides a wide-band, diffusive radiation
    condition, where the patch temperature is specified.

Usage
    \table
        Property     | Description             | Required    | Default value
        T            | temperature field name  | no          | T
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            wideBandDiffusiveRadiation;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField
    Foam::radiationCoupledBase

SourceFiles
    wideBandDiffusiveRadiationMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wideBandDiffusiveRadiationMixedFvPatchScalarField_H
#define wideBandDiffusiveRadiationMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "radiationCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class wideBandDiffusiveRadiationMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wideBandDiffusiveRadiationMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public radiationCoupledBase
{
    // Private Data

        //- Name of temperature field
        word TName_;


public:

    //- Runtime type information
    TypeName("wideBandDiffusiveRadiation");


    // Constructors

        //- Construct from patch and internal field
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given GreyDiffusiveRadiationMixedFvPatchField
        //  onto a new patch
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wideBandDiffusiveRadiationMixedFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationMixedFvPatchScalarField&,
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
                new wideBandDiffusiveRadiationMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the temperature field name
            const word& TName() const
            {
                return TName_;
            }

            //- Return reference to the temperature field name to allow
            //  adjustment
            word& TName()
            {
                return TName_;
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


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

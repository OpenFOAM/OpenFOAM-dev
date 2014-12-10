/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a grey-diffuse condition for radiation
    intensity, \c I, for use with the finite-volume discrete-ordinates model
    (fvDOM), in which the radiation temperature is retrieved from the
    temperature field boundary condition.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        T            | temperature field name  | no          | T
        emissivityMode | emissivity mode: solidThermo or lookup | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            greyDiffusiveRadiation;
        T               T;
        emissivityMode  solidThermo;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::radiation::radiationModel
    Foam::radiation::fvDOM
    Foam::radiationCoupledBase
    Foam::mixedFvPatchField

SourceFiles
    greyDiffusiveRadiationMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef greyDiffusiveRadiationMixedFvPatchScalarField_H
#define greyDiffusiveRadiationMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "radiationCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
        Class greyDiffusiveRadiationMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class greyDiffusiveRadiationMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public radiationCoupledBase
{
    // Private data

        //- Name of temperature field
        word TName_;


public:

    //- Runtime type information
    TypeName("greyDiffusiveRadiation");


    // Constructors

        //- Construct from patch and internal field
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  greyDiffusiveRadiationMixedFvPatchScalarField onto a new patch
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new greyDiffusiveRadiationMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationMixedFvPatchScalarField&,
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
                new greyDiffusiveRadiationMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

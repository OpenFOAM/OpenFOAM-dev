/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::adsorptionMassFractionFvPatchScalarField

Description
    This is a mass-fraction boundary condition for an adsorbing wall.

    This condition models a baffle which adsorbs some species and not others.
    It must be used in conjunction with a
    specieTransferVelocityFvPatchVectorField velocity condition, and a
    specieTransferTemperatureFvPatchScalarField  temperature condition.

    The mass flux of a species is calculated as a coefficient multiplied by the
    concentration adjacent to the patch.
        \f[
            \phi_{Yi} = c A \psi_i
        \f]
        where
        \vartable
            \phi_{Yi}  | Flux of the adsorbing specie [kg/s]
            c          | Transfer coefficient [kg/m^2/s/<property-dimensions>]
            A          | Patch face area [m^2]
            C_i        | Property adjacent to the patch [<property-dimensions>]
        \endvartable

    A species that is adsorbed by the patch will, therefore, have a coefficient
    greater than zero, whilst a species that is not adsorbed will have a
    coefficient equal to zero.

Usage
    \table
        Property | Description                | Req'd? | Default
        c        | Transfer coefficient       | no     | 0
        property | Property used to drive the transfer; massFraction, \\
        moleFraction, molarConcentration, or partialPressure | if c is \\
        non-zero |
        phi      | Name of the flux field     | no     | phi
        U        | Name of the velocity field | no     | U
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            adsorptionMassFraction;
        property        molarConcentration;
        c               1e-3;
        value           $internalField;
    }
    \endverbatim

See also
    Foam::specieTransferMassFractionFvPatchScalarField
    Foam::specieTransferVelocityFvPatchVectorField
    Foam::specieTransferTemperatureFvPatchScalarField

SourceFiles
    adsorptionMassFractionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adsorptionMassFractionFvPatchScalarField_H
#define adsorptionMassFractionFvPatchScalarField_H

#include "specieTransferMassFractionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicSpecieMixture;

/*---------------------------------------------------------------------------*\
          Class adsorptionMassFractionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class adsorptionMassFractionFvPatchScalarField
:
    public specieTransferMassFractionFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("adsorptionMassFraction");


    // Constructors

        //- Construct from patch and internal field
        adsorptionMassFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adsorptionMassFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        adsorptionMassFractionFvPatchScalarField
        (
            const adsorptionMassFractionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Disallow copy without setting internal field reference
        adsorptionMassFractionFvPatchScalarField
        (
            const adsorptionMassFractionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        adsorptionMassFractionFvPatchScalarField
        (
            const adsorptionMassFractionFvPatchScalarField&,
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
                new adsorptionMassFractionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Evaluation functions

            //- Return the flux of this species
            tmp<scalarField> calcPhiYp() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Foam::semiPermeableBaffleMassFractionFvPatchScalarField

Description
    This is a mass-fraction boundary condition for a semi-permeable baffle.

    This condition models a baffle which is permeable to a some species and
    impermeable to others. It must be used in conjunction with a
    specieTransferVelocityFvPatchVectorField velocity condition, and a
    specieTransferTemperatureFvPatchScalarField  temperature condition.

    The mass flux of a species is calculated as a coefficient multiplied by the
    difference in a property across the baffle.
        \f[
            \phi_{Yi} = c A (\psi_i - \psi_{i,n})
        \f]
        where
        \vartable
            \phi_{Yi}  | Flux of the permeable specie [kg/s]
            c          | Transfer coefficient [kg/m^2/s/<property-dimensions>]
            A          | Patch face area [m^2]
            \psi_i     | Property on the patch [<property-dimensions>]
            \psi_{i,n} | Property on the neighbour patch [<property-dimensions>]
        \endvartable

    A species that the baffle is permeable to will, therefore, have a
    coefficient greater than zero, whilst a species that does not transfer will
    have a coefficient equal to zero.

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
        type            semiPermeableBaffleMassFraction;
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
    semiPermeableBaffleMassFractionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef semiPermeableBaffleMassFractionFvPatchScalarField_H
#define semiPermeableBaffleMassFractionFvPatchScalarField_H

#include "specieTransferMassFractionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class semiPermeableBaffleMassFractionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class semiPermeableBaffleMassFractionFvPatchScalarField
:
    public specieTransferMassFractionFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("semiPermeableBaffleMassFraction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        semiPermeableBaffleMassFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        semiPermeableBaffleMassFractionFvPatchScalarField
        (
            const semiPermeableBaffleMassFractionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        semiPermeableBaffleMassFractionFvPatchScalarField
        (
            const semiPermeableBaffleMassFractionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        semiPermeableBaffleMassFractionFvPatchScalarField
        (
            const semiPermeableBaffleMassFractionFvPatchScalarField&,
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
                new semiPermeableBaffleMassFractionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Evaluation functions

            //- Return the flux of this species through the baffle
            virtual tmp<scalarField> calcPhiYp() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

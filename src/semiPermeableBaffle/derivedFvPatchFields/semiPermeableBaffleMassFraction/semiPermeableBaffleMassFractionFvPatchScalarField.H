/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    impermeable to others. It must be used in conjunction with the
    corresponding velocity condition,
    semiPermeableBaffleVelocityFvPatchVectorField.

    The mass flux of a species is calculated as a coefficient multiplied by the
    difference in an input variable across the baffle.
        \f[
            \phi_{Yi} = c A (\psi_i - \psi_{i,n})
        \f]
        where
        \vartable
            \phi_{Yi}  | flux of the permeable species [kg/s]
            c          | transfer coefficient [kg/m^2/s/<input-dimensions>]
            A          | patch face area [m^2]
            \psi_i     | input variable on the patch []
            \psi_{i,n} | input variable on the neighbour patch []
        \endvartable

    A species that the baffle is permeable to will, therefore, have a
    coefficient greater than zero, whilst a species that does not transfer will
    have a coefficient equal to zero.

    This condition calculates the species flux. The fluxes are summed up by the
    velocity condition to generate the net mass transfer across the baffle.
    This mass-fraction condition then generates a corrective diffusive flux to
    ensure that the correct amounts of the permeable species are transferred.

Usage
    \table
        Property | Description                | Req'd? | Default
        c        | Transfer coefficient       | no     | 0
        input    | Input variable used to drive the transfer; massFraction, \\
        moleFraction or partialPressure | if c is non-zero | none
        phi      | Name of the flux field     | no     | phi
        p        | Name of the pressure field | no     | p
    \endtable

See also
    Foam::semiPermeableBaffleVelocityFvPatchVectorField

SourceFiles
    semiPermeableBaffleMassFractionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef semiPermeableBaffleMassFractionFvPatchScalarField_H
#define semiPermeableBaffleMassFractionFvPatchScalarField_H

#include "mappedPatchBase.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicSpecieMixture;

/*---------------------------------------------------------------------------*\
      Class semiPermeableBaffleMassFractionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class semiPermeableBaffleMassFractionFvPatchScalarField
:
    public mappedPatchBase,
    public mixedFvPatchScalarField
{
public:

        //- Enumeration for the input variable driving the transfer
        enum input
        {
            none,
            massFraction,
            moleFraction,
            partialPressure,
        };

        //- Input variable type names
        static const NamedEnum<input, 4> inputNames_;


private:

    // Private Data

        //- Transfer coefficient
        const scalar c_;

        //- Input variable driving the transfer
        const input input_;

        //- Name of the flux field
        const word phiName_;

        //- Name of the pressure field. Only needed if mode is pressure.
        const word pName_;


public:

    //- Runtime type information
    TypeName("semiPermeableBaffleMassFraction");


    // Static member functions

        //- Access the composition for the given database
        static const basicSpecieMixture& composition(const objectRegistry& db);


    // Constructors

        //- Construct from patch and internal field
        semiPermeableBaffleMassFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

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
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        semiPermeableBaffleMassFractionFvPatchScalarField
        (
            const semiPermeableBaffleMassFractionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new semiPermeableBaffleMassFractionFvPatchScalarField(*this)
            );
        }

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
            tmp<scalarField> phiY() const;

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

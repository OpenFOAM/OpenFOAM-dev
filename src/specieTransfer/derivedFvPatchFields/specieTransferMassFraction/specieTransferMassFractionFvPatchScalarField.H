/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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
    Foam::specieTransferMassFractionFvPatchScalarField

Description
    Abstract base class for specie-transferring mass fraction boundary
    conditions.

    Derived classes compute the actual species flux. Those fluxes are summed up
    by the specieTransferVelocityFvPatchVectorField velocity condition to
    generate the net mass transfer into or out of the patch. This mass-fraction
    condition then generates a corrective diffusive flux to ensure that the
    correct amounts of the transferring species are transported through the
    patch.

SourceFiles
    specieTransferMassFractionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef specieTransferMassFractionFvPatchScalarField_H
#define specieTransferMassFractionFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidMulticomponentThermo;

/*---------------------------------------------------------------------------*\
        Class specieTransferMassFractionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class specieTransferMassFractionFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

        //- Enumeration for the property driving the transfer
        enum property
        {
            massFraction,
            moleFraction,
            molarConcentration,
            partialPressure
        };

        //- Property type names
        static const NamedEnum<property, 4> propertyNames_;


private:

    // Private Data

        //- Name of the flux field
        const word phiName_;

        //- Name of the velocity field
        const word UName_;

        //- Flux of this species
        mutable scalarField phiYp_;

        //- Current time index
        mutable label timeIndex_;


protected:

    // Protected Data

        //- Transfer coefficient
        const scalar c_;

        //- Property driving the transfer
        const property property_;


public:

    // Static member functions

        //- Get the thermo from the given database
        static const fluidMulticomponentThermo& thermo
        (
            const objectRegistry& db
        );


    // Constructors

        //- Construct from patch, internal field and dictionary
        specieTransferMassFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  specieTransferMassFractionFvPatchScalarField
        //  onto a new patch
        specieTransferMassFractionFvPatchScalarField
        (
            const specieTransferMassFractionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        specieTransferMassFractionFvPatchScalarField
        (
            const specieTransferMassFractionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        specieTransferMassFractionFvPatchScalarField
        (
            const specieTransferMassFractionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Return the flux of this species
            virtual tmp<scalarField> calcPhiYp() const = 0;

            //- Return the flux of this species
            virtual const scalarField& phiYp() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


    // Member Operators

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using mixedFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::lumpedMassTemperatureFvPatchScalarField

Description
    This boundary condition is applied to a patch which bounds a solid body,
    wholly or partially.  It represents the body as a lumped mass, i.e. by a
    single temperature \c T which is fixed across the patch.  The body has a
    volume \c V which is either specified by the user, or is calculated when
    the patch describes a closed volume (including in 2D meshes). Starting from
    an initial, specified \c T, the change in temperature is calculated over
    time according to an applied power source \c Q and the heat transferred
    across the boundary \f$Q_{b}\f$ (positive into the lumped mass):

        \f[
             dT/dt = frac{Q + Q_{b}}{\rho C_{v} V}
        \f]

    where
    \vartable
        Q      | specified power source [W]
        Q_{b}  | total calculated heat transferred across the boundary [W]
        \rho   | density [kg/m^3]
        C_{v}  | specific heat capacity [J/(kg K)]
        V      | volume of the lumped mass [m^3}
    \endtable

Usage
    \table
        Property   | Description                  | Req'd? | Default
        rho        | density                      | yes    |
        Cv         | specific heat capacity       | yes    |
        T          | temperature                  | yes    |
        Q          | power source                 | no     | 0
        volume     | volume of the lumped mass    | no     | calculated
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       lumpedMassTemperature;
        rho        1000;
        Cv         1300;
        T          500;
        Q          constant 0.5;
        value      uniform $T; // placeholder
    }
    \endverbatim

SourceFiles
    lumpedMassTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedMassTemperatureFvPatchScalarField_H
#define lumpedMassTemperatureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "UniformDimensionedField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class lumpedMassTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class lumpedMassTemperatureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Density
        const scalar rho_;

        //- Specific heat capacity
        const scalar Cv_;

        //- Temperature of the thermal mass
        UniformDimensionedField<scalar> T_;

        //- Heat source (optional)
        autoPtr<Function1<scalar>> Q_;

        //- Volume of the lumped mass, defaults when enclosed by the patch
        scalar V_;


    // Private Member Functions

        //- Does the patch form a closed volume, including in a 2D case?
        bool closed() const;

        //- Return volume formed by the enclosed patch
        scalar V() const;


public:

    //- Runtime type information
    TypeName("lumpedMassTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        lumpedMassTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        lumpedMassTemperatureFvPatchScalarField
        (
            const lumpedMassTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        lumpedMassTemperatureFvPatchScalarField
        (
            const lumpedMassTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        lumpedMassTemperatureFvPatchScalarField
        (
            const lumpedMassTemperatureFvPatchScalarField&,
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
                new lumpedMassTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


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

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
    Foam::PrghPressureFvPatchScalarField

Description
    This boundary condition provides the p_rgh equivalent of a pressure
    boundary condition calculated as:

        \f[
            p_rgh = p - \rho g (h - hRef)
        \f]

    where
    \vartable
        p_rgh   | Pseudo hydrostatic pressure [Pa]
        p       | Static pressure [Pa]
        h       | Height in the opposite direction to gravity
        hRef    | Reference height in the opposite direction to gravity
        \rho    | density
        g       | acceleration due to gravity [m/s^2]
    \endtable

SourceFiles
    PrghPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef PrghPressureFvPatchScalarField_H
#define PrghPressureFvPatchScalarField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class PrghPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

template<class PressureFvPatchScalarField>
class PrghPressureFvPatchScalarField
:
    public PressureFvPatchScalarField
{
    // Private Data

        //- Name of the density field, defaults to "rho"
        const word rhoName_;


public:

    //- Runtime type information
    TypeName("PrghPressure");


    // Constructors

        //- Construct from patch and internal field
        PrghPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        PrghPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  PrghPressureFvPatchScalarField onto a new patch
        PrghPressureFvPatchScalarField
        (
            const PrghPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        PrghPressureFvPatchScalarField
        (
            const PrghPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        PrghPressureFvPatchScalarField
        (
            const PrghPressureFvPatchScalarField&,
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
                new PrghPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

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
    #include "PrghPressureFvPatchScalarField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePrghPatchScalarField(Pressure, PrghPressure)                       \
    typedef PrghPressureFvPatchScalarField<Pressure##FvPatchScalarField>       \
    PrghPressure;                                                              \
                                                                               \
    defineTemplateTypeNameAndDebug(PrghPressure, 0);                           \
                                                                               \
    addToPatchFieldRunTimeSelection                                            \
    (                                                                          \
        fvPatchScalarField,                                                    \
        PrghPressure                                                           \
    )

#endif

// ************************************************************************* //

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
    Foam::prghCyclicPressureFvPatchScalarField

Description
    This boundary condition provides a cyclic condition for p_rgh. It applies
    corrections to the value and gradient on both sides of the cyclic to
    account for the non-cylicity of the gravitational force.

    This condition is only needed when the cyclic patches have a transformation
    and a normal component in the direction of gravity. If the cyclic patches
    are orthogonal to the direction gravity, then a normal cyclic boundary
    condition can be used instead.

    Care must be taken when using this boundary condition that the simulation
    is actually cyclic. The following constraints apply:

    - Both cyclic patches must be oriented in the same way with respect to
      gravity. In practice this means that applicability is limited to cyclics
      with translational transformations.

    - The model cannot have any dependence on the absolute value of the
      pressure field. The absolute value of the pressure, in reality, varies
      between each repetition of the geometry; it is not actually formally
      cyclic. Only the gradient of the pressure field can be truly cyclic. This
      model is therefore only valid if the absolute value of the pressure is
      arbitrary, and only the gradient has an effect on the solution. This is
      the case for incompressible multiphase solutions or incompressible
      Boussinesq-like models of density variation. It is not true if (for
      example) a compressible thermodynamic model is being used.

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type (should be \c cyclic) | yes |
        rhoInf       | far-field density       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            prghCyclicPressure;
        patchType       cyclic;
        rhoInf          1;
    }
    \endverbatim

SourceFiles
    prghCyclicPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef prghCyclicPressureFvPatchScalarField_H
#define prghCyclicPressureFvPatchScalarField_H

#include "jumpCyclicFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class prghCyclicPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class prghCyclicPressureFvPatchScalarField
:
    public jumpCyclicFvPatchScalarField
{
    // Private Data

        //- Name of the density field
        const word rhoName_;

        //- Far-field density
        const scalar rhoInf_;

        //- Jump in value from the other patch to this one
        scalarField jump_;


    // Private Member Functions

        //- Get the far-field density from the owner patch
        scalar rhoInf() const;

        //- Update result field based on interface functionality
        template<class ... Cmpt>
        void updateInterfaceMatrix
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const Pstream::commsTypes commsType,
            const Cmpt ... cmpt
        ) const;


public:

    //- Runtime type information
    TypeName("prghCyclicPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        prghCyclicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        prghCyclicPressureFvPatchScalarField
        (
            const prghCyclicPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        prghCyclicPressureFvPatchScalarField
        (
            const prghCyclicPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        prghCyclicPressureFvPatchScalarField
        (
            const prghCyclicPressureFvPatchScalarField&,
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
                new prghCyclicPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return the "jump"
            virtual tmp<scalarField> jump() const;


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the patch pressure gradient field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

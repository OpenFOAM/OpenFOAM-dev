/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::porousBafflePressureFvPatchField

Description

    This boundary condition provides a porous baffle pressure jump condition,
    using the \c cyclic condition as a base. The jump in pressure is defined
    by:

        \f[
            \Delta p = -(D \mu U + 0.5 I \rho |U|^2 )L
        \f]

    where

    \vartable
        p      | pressure [Pa]
        \rho   | density [kg/m^3]
        \mu    | laminar viscosity [Pa s]
        D      | Darcy coefficient
        I      | inertial coefficient
        L      | porous media length in the flow direction
    \endvartable


Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type (should be \c cyclic) | yes |
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        D            | Darcy coefficient       | yes         |
        I            | inertial coefficient    | yes         |
        length       | porous media length in the flow direction | yes |
        relaxation   | relaxation factor for pressure jump | no | 1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            porousBafflePressure;
        patchType       cyclic;
        jump            uniform 0;
        D               0.001;
        I               1000000;
        length          0.1;
        relaxation      0.2;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    porousBafflePressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef porousBafflePressureFvPatchField_H
#define porousBafflePressureFvPatchField_H

#include "fixedJumpFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class porousBafflePressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class porousBafflePressureFvPatchField
:
    public fixedJumpFvPatchScalarField
{
    // Private Data

        //- Name of flux field (default = phi)
        const word phiName_;

        //- Name of density field (default = rho)
        const word rhoName_;

        //- Darcy pressure loss coefficient
        const scalar D_;

        //- Inertia pressure lost coefficient
        const scalar I_;

        //- Porous media length
        const scalar length_;

        //- Relaxation for pressure jump
        const scalar relaxation_;

        //- Previous iteration "jump" field for relaxation
        scalarField jump0_;


public:

    //- Runtime type information
    TypeName("porousBafflePressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        porousBafflePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given porousBafflePressureFvPatchField
        // onto a new patch
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField&
        ) = delete;

        //- Copy constructor setting internal field reference
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new porousBafflePressureFvPatchField(*this, iF)
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

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::fluxCorrectedVelocityFvPatchVectorField

Description
    This boundary condition provides a velocity outlet boundary condition for
    patches where the pressure is specified.  The outflow velocity is obtained
    by "zeroGradient" and then corrected from the flux:

        \f[
            U_p = U_c - n (n \cdot U_c) + \frac{n \phi_p}{|S_f|}
        \f]

    where

    \vartable
        U_p | velocity at the patch [m/s]
        U_c | velocity in cells adjacent to the patch [m/s]
        n   | patch normal vectors
        \phi_p | flux at the patch [m^3/s or kg/s]
        S_f | patch face area vectors [m^2]
    \endvartable

    where

    \table
        Property     | Description             | Required    | Default value
        phi          | name of flux field      | no          | phi
        rho          | name of density field   | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fluxCorrectedVelocity;
        phi             phi;
        rho             rho;
    }
    \endverbatim

    Note:
        If reverse flow is possible or expected use the
        pressureInletOutletVelocity condition instead.

See also
    Foam::zeroGradientFvPatchField
    Foam::pressureInletOutletVelocityFvPatchVectorField

SourceFiles
    fluxCorrectedVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fluxCorrectedVelocityFvPatchVectorField_H
#define fluxCorrectedVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class fluxCorrectedVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fluxCorrectedVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // Private Data

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;


public:

    //- Runtime type information
    TypeName("fluxCorrectedVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        fluxCorrectedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fluxCorrectedVelocityFvPatchVectorField
        //  onto a new patch
        fluxCorrectedVelocityFvPatchVectorField
        (
            const fluxCorrectedVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        fluxCorrectedVelocityFvPatchVectorField
        (
            const fluxCorrectedVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new fluxCorrectedVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

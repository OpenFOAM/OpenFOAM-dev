/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::maxwellSlipUFvPatchVectorField

Description
    Maxwell slip boundary condition including thermal creep and surface
    curvature terms that can be optionally switched off.

SourceFiles
    fixedRhoFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef maxwellSlipUFvPatchVectorField_H
#define maxwellSlipUFvPatchVectorField_H

#include "mixedFixedValueSlipFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class maxwellSlipUFvPatch Declaration
\*---------------------------------------------------------------------------*/

class maxwellSlipUFvPatchVectorField
:
    public mixedFixedValueSlipFvPatchVectorField
{
    // Private data

        //- Temperature field name, default = "T"
        word TName_;

        //- Density field name, default = "rho"
        word rhoName_;

        //- Compressibility field name, default = "thermo:psi"
        word psiName_;

        //- Dynamic viscosity field name, default = "thermo:mu"
        word muName_;

        //- tauMC field name, default = "tauMC"
        word tauMCName_;

        // Accommodation coefficient
        scalar accommodationCoeff_;

        // Wall velocity
        vectorField Uwall_;

        // Include thermal creep term (default on)
        Switch thermalCreep_;

        // Include boundary curvature term (default on)
        Switch curvature_;


public:

    //- Runtime type information
    TypeName("maxwellSlipU");


    // Constructors

        //- Construct from patch and internal field
        maxwellSlipUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        maxwellSlipUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  maxwellSlipUFvPatchVectorField onto a new patch
        maxwellSlipUFvPatchVectorField
        (
            const maxwellSlipUFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new maxwellSlipUFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        maxwellSlipUFvPatchVectorField
        (
            const maxwellSlipUFvPatchVectorField&,
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
                new maxwellSlipUFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::gradientEnergyFvPatchScalarField

Description
    This boundary condition provides a gradient condition for internal energy,
    where the gradient is calculated using:

        \f[
            \nabla(e_p) = \nabla_\perp C_p(p, T) + \frac{e_p - e_c}{\Delta}
        \f]

    where
    \vartable
        e_p     | energy at patch faces [J]
        e_c     | energy at patch internal cells [J]
        p       | pressure [bar]
        T       | temperature [K]
        C_p     | specific heat [J/kg/K]
        \Delta  | distance between patch face and internal cell centres [m]
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            gradientEnergy;
        gradient        uniform 10;
    }
    \endverbatim

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    gradientEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef gradientEnergyFvPatchScalarField_H
#define gradientEnergyFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class gradientEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class gradientEnergyFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("gradientEnergy");


    // Constructors

        //- Construct from patch and internal field
        gradientEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gradientEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given gradientEnergyFvPatchScalarField
        // onto a new patch
        gradientEnergyFvPatchScalarField
        (
            const gradientEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        gradientEnergyFvPatchScalarField
        (
            const gradientEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new gradientEnergyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gradientEnergyFvPatchScalarField
        (
            const gradientEnergyFvPatchScalarField&,
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
                new gradientEnergyFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

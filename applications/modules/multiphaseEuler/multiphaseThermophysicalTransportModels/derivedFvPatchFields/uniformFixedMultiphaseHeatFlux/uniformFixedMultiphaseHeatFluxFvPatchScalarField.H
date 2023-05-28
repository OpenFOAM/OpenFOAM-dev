/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
    Foam::uniformFixedMultiphaseHeatFluxFvPatchScalarField

Description
    Uniform fixed heat flux boundary condition for Eulerian multi-phase cases.
    Constructs a mixed constraint which portions the heat flux between the
    phases in such a way as to keep the boundary temperature uniform across all
    phases. The heat flux can be specified as a time-varying function, and an
    under-relaxation factor can be supplied if this is necessary to maintain
    stability.

Usage
    \table
        Property     | Description                | Required | Default value
        q            | Heat flux [w/m^2]          | yes      |
        relax        | Relaxation factor          | no       | 1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformFixedMultiphaseHeatFlux;
        q               1000;
        relax           0.3;
        value           $internalField;
    }
    \endverbatim

SourceFiles
    uniformFixedMultiphaseHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformFixedMultiphaseHeatFluxFvPatchScalarField_H
#define uniformFixedMultiphaseHeatFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class uniformFixedMultiphaseHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class uniformFixedMultiphaseHeatFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Heat flux [W/m^2]
        autoPtr<Function1<scalar>> q_;

        //- Relaxation factor
        scalar relax_;


public:

    //- Runtime type information
    TypeName("uniformFixedMultiphaseHeatFlux");


    // Constructors

        //- Construct from patch, internal field and dictionary
        uniformFixedMultiphaseHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  uniformFixedMultiphaseHeatFluxFvPatchScalarField
        //  onto a new patch
        uniformFixedMultiphaseHeatFluxFvPatchScalarField
        (
            const uniformFixedMultiphaseHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        uniformFixedMultiphaseHeatFluxFvPatchScalarField
        (
            const uniformFixedMultiphaseHeatFluxFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        uniformFixedMultiphaseHeatFluxFvPatchScalarField
        (
            const uniformFixedMultiphaseHeatFluxFvPatchScalarField&,
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
                new uniformFixedMultiphaseHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::compressible::alphatJayatillekeWallFunctionFvPatchScalarField

Description
    This boundary condition provides a thermal wall function for turbulent
    thermal diffusivity (usually\c alphat) based on the Jayatilleke model.

Usage
    \table
        Property     | Description             | Required    | Default value
        Prt          | turbulent Prandtl number | no         | 0.85
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            alphatJayatillekeWallFunction;
        Prt             0.85;
        value           uniform 0;
    }
    \endverbatim

    Note that other model constants (i.e., Cmu, kappa and E) are obtained from
    the corresponding turbulent viscosity boundary condition.

See also
    Foam::fixedValueFvPatchField

SourceFiles
    alphatJayatillekeWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatJayatillekeWallFunctionFvPatchScalarField_H
#define alphatJayatillekeWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidThermophysicalTransportModel;

namespace compressible
{

/*---------------------------------------------------------------------------*\
       Class alphatJayatillekeWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatJayatillekeWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Static Data

        //- Solution tolerance
        static const scalar tolerance_;

        //- Maximum number of solution iterations
        static const label maxIters_;


    // Private Data

        //- Turbulent Prandtl number
        scalar Prt_;


public:

    //- Runtime type information
    TypeName("compressible::alphatJayatillekeWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given an
        //  alphatJayatillekeWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionFvPatchScalarField&,
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
                new alphatJayatillekeWallFunctionFvPatchScalarField
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


        // Evaluation functions

            //- Calculate the smoothing function
            static tmp<scalarField> P(const scalarField& Prat);

            //- Calculate y+ at the edge of the thermal laminar sublayer
            static tmp<scalarField> yPlusTherm
            (
                const nutWallFunctionFvPatchScalarField& nutw,
                const scalarField& P,
                const scalarField& Prat
            );

            //- Calculate the turbulent thermal diffusivity
            static tmp<scalarField> alphat
            (
                const fluidThermophysicalTransportModel& ttm,
                const scalar Prt,
                const label patchi
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::omegaWallFunctionFvPatchScalarField

Description
    This boundary condition provides a wall constraint on turbulnce specific
    dissipation, omega for both low and high Reynolds number turbulence models.

    The near-wall omega may be either blended between the viscous region and
    logarithmic region values or switched between these values based on the
    laminar-to-turbulent y+ value derived from kappa and E specified in the
    corresponding nutWallFunction.  Recent tests have shown that the standard
    switching method provides more accurate results for 10 < y+ < 30 when used
    with high Reynolds number wall-functions and both methods provide accurate
    results when used with continuous wall-functions.  Based on this the
    standard switching method is used by default.

Usage
    \table
        Property     | Description             | Required    | Default value
        beta1        | Model coefficient       | no          | 0.075
        blended      | Blending switch         | no          | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            omegaWallFunction;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField
    Foam::epsilonWallFunctionFvPatchScalarField

SourceFiles
    omegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaWallFunctionFvPatchScalarField_H
#define omegaWallFunctionFvPatchScalarField_H

#include "wallCellWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class momentumTransportModel;

/*---------------------------------------------------------------------------*\
             Class omegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaWallFunctionFvPatchScalarField
:
    public wallCellWallFunctionFvPatchScalarField
{
private:

    // Private Data

        //- beta1 coefficient
        scalar beta1_;

        //- Blending switch (defaults to false)
        Switch blended_;

        //- Wall cell turbulence generation values. Master patch only.
        autoPtr<scalarField> wallCellGPtr_;

        //- Wall cell turbulence frequency values. Master patch only.
        autoPtr<scalarField> wallCellOmegaPtr_;


    // Private Member Functions

        //- Calculate and return G and omega
        Pair<tmp<scalarField>> calculate(const momentumTransportModel&) const;

        //- Calculate the turbulent values in wall cells and apply to the
        //  turbulence model's fields
        void updateCoeffsMaster();

        //- Constrain the matrix such that the wall cells take the value
        //  specified by the wall function
        void manipulateMatrixMaster(fvMatrix<scalar>& matrix);


public:

    //- Runtime type information
    TypeName("omegaWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        omegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  omegaWallFunctionFvPatchScalarField
        //  onto a new patch
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
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
                new omegaWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~omegaWallFunctionFvPatchScalarField()
    {}


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

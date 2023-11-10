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
    Foam::nutWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on turbulence kinetic energy.
    - replicates OpenFOAM v1.5 (and earlier) behaviour

Usage
    \table
        Property  | Description         | Required   | Default value
        Cmu       | Cmu coefficient     | no         | 0.09
        kappa     | Von Karman constant | no         | 0.41
        E         | E coefficient       | no         | 9.8
    \endtable

    Examples of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutWallFunction;
        value           uniform 0;
    }
    \endverbatim

    Reference for the default model coefficients:
    \verbatim
        H. Versteeg, W. Malalasekera
        An Introduction to Computational Fluid Dynamics: The Finite Volume
        Method, subsection "3.5.2 k-epsilon model"
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    nutWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutWallFunctionFvPatchScalarField_H
#define nutWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class momentumTransportModel;

/*---------------------------------------------------------------------------*\
              Class nutWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const = 0;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        nutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutWallFunctionFvPatchScalarField
        //  onto a new patch
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions

        //- Return Cmu
        scalar Cmu() const
        {
            return Cmu_;
        }

        //- Return kappa
        scalar kappa() const
        {
            return kappa_;
        }

        //- Return E
        scalar E() const
        {
            return E_;
        }

        //- Return the nut patchField for the given wall patch
        static const nutWallFunctionFvPatchScalarField& nutw
        (
            const momentumTransportModel& turbModel,
            const label patchi
        );

        //- Calculate the Y+ at the edge of the laminar sublayer
        static scalar yPlusLam(const scalar kappa, const scalar E);

        //- Return the Y+ at the edge of the laminar sublayer
        scalar yPlusLam() const;

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const = 0;


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

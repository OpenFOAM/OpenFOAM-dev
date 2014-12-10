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
    Foam::compressible::mutWallFunctionFvPatchScalarField

Group
    grpCmpWallFunctions

Description
    This boundary condition provides a turbulent viscosity condition when
    using wall functions, based on turbulence kinetic energy.
    - replicates OpenFOAM v1.5 (and earlier) behaviour

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            mutWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    mutWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutWallFunctionFvPatchScalarField_H
#define compressibleMutWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
             Class mutWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutWallFunctionFvPatchScalarField
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
        virtual tmp<scalarField> calcMut() const = 0;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("mutWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutWallFunctionFvPatchScalarField
        //  onto a new patch
        mutWallFunctionFvPatchScalarField
        (
            const mutWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutWallFunctionFvPatchScalarField
        (
            const mutWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        mutWallFunctionFvPatchScalarField
        (
            const mutWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions

        //- Calculate the Y+ at the edge of the laminar sublayer
        static scalar yPlusLam(const scalar kappa, const scalar E);

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

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

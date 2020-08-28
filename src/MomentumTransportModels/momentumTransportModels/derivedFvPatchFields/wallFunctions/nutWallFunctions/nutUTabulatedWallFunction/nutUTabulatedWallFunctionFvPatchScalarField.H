/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::nutUTabulatedWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions.  As input, the user specifies a look-up table
    of U+ as a function of near-wall Reynolds number.  The table should be
    located in the $FOAM_CASE/constant directory.

Usage
    \table
        Property     | Description             | Required    | Default value
        uPlusTable   | U+ as a function of Re table name | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutTabulatedWallFunction;
        uPlusTable      myUPlusTable;
        value           uniform 0;
    }
    \endverbatim

Note
    The tables are not registered since the same table object may be used for
    more than one patch.

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutUTabulatedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUTabulatedWallFunctionFvPatchScalarField_H
#define nutUTabulatedWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"
#include "uniformInterpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class nutUTabulatedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUTabulatedWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- Name of u+ table
        word uPlusTableName_;

        //- U+ table
        uniformInterpolationTable<scalar> uPlusTable_;


    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;

        //- Calculate wall u+ from table
        virtual tmp<scalarField> calcUPlus(const scalarField& Rey) const;


public:

    //- Runtime type information
    TypeName("nutTabulatedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUTabulatedWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUTabulatedWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&,
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
                new nutUTabulatedWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

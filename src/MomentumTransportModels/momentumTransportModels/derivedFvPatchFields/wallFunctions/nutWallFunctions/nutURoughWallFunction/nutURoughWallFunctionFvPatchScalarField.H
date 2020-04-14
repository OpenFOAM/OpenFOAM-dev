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
    Foam::nutURoughWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions for rough walls, based on velocity.  The
    condition manipulates the E parameter to account for roughness effects.

    Reference:
    \verbatim
        Cebeci, T., & Bradshaw, P. (1977).
        Momentum transfer in boundary layers.
        Washington, DC, Hemisphere Publishing Corp.; New York,
        McGraw-Hill Book Co., 1977. 407 p.
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        Ks           | sand-grain roughness height | yes     |
        Cs           | roughness constant      | yes         |
    \endtable

    Parameter ranges
    - roughness height (Ks) = sand-grain roughness (0 for smooth walls)
    - roughness constant (Cs) = 0.5 - 1.0

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutURoughWallFunction;
        Ks              uniform 1e-5;
        Cs              uniform 0.5;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::nutUWallFunctionFvPatchScalarField

SourceFiles
    nutURoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutURoughWallFunctionFvPatchScalarField_H
#define nutURoughWallFunctionFvPatchScalarField_H

#include "nutUWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutURoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutURoughWallFunctionFvPatchScalarField
:
    public nutUWallFunctionFvPatchScalarField
{
    // Private Data

        //- Roughness height
        scalarField Ks_;

        //- Roughness constant
        scalarField Cs_;


protected:

    // Protected Member Functions

        //- Calculate yPLus
        virtual tmp<scalarField> yPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("nutURoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutURoughWallFunctionFvPatchScalarField
        //  onto a new patch
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutURoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&,
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
                new nutURoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access functions

            //- Return the roughness height
            const scalarField& Ks() const
            {
                return Ks_;
            }

            //- Return reference to the roughness height to allow adjustment
            scalarField& Ks()
            {
                return Ks_;
            }


            //- Return the roughness constant scale
            const scalarField& Cs() const
            {
                return Cs_;
            }

            //- Return reference to the roughness constant to allow adjustment
            scalarField& Cs()
            {
                return Cs_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

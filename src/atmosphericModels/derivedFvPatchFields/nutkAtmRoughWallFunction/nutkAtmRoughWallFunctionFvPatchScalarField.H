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
    Foam::nutkAtmRoughWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent kinematic viscosity for
    atmospheric velocity profiles.  It is designed to be used in conjunction
    with the atmBoundaryLayerInletVelocity boundary condition.  The values
    are calculated using:

        \f[
            U = frac{U_f}{\kappa} ln(\frac{z + z_0}{z_0})
        \f]

    where

    \vartable
        U_f    | frictional velocity
        \kappa | Von Karman's constant
        z_0    | surface roughness length
        z      | vertical co-ordinate
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        z0           | surface roughness length| yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutkAtmRoughWallFunction;
        z0              uniform 0;
    }
    \endverbatim

See also
    Foam::nutkWallFunctionFvPatchField

SourceFiles
    nutkAtmRoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkAtmRoughWallFunctionFvPatchScalarField_H
#define nutkAtmRoughWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class nutkAtmRoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkAtmRoughWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- Surface roughness length
        scalarField z0_;


    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("nutkAtmRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkAtmRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkAtmRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutkAtmRoughWallFunctionFvPatchScalarField
        //  onto a new patch
        nutkAtmRoughWallFunctionFvPatchScalarField
        (
            const nutkAtmRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        nutkAtmRoughWallFunctionFvPatchScalarField
        (
            const nutkAtmRoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutkAtmRoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        nutkAtmRoughWallFunctionFvPatchScalarField
        (
            const nutkAtmRoughWallFunctionFvPatchScalarField&,
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
                new nutkAtmRoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access functions

            // Return z0
            scalarField& z0()
            {
                return z0_;
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
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

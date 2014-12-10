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
    Foam::compressible::RASModels::alphatFilmWallFunctionFvPatchScalarField

Group
    grpSurfaceFilmBoundaryConditions grpCmpWallFunctions

Description
    This boundary condition provides a turbulent thermal diffusivity condition
    when using wall functions, for use with surface film models.  This
    condition varies from the standard wall function by taking into account any
    mass released from the film model.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        B            | model coefficient       | no          | 5.5
        yPlusCrit    | critical y+ for transition to turbulent flow | no|11.05
        Cmu          | model coefficient       | no          | 0.09
        kappa        | Von-Karman constant     | no          | 0.41
        Prt          | turbulent Prandtl number | no         | 0.85
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            alphatFilmWallFunction;
        B               5.5;
        yPlusCrit       11.05;
        Cmu             0.09;
        kappa           0.41;
        Prt             0.85;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    alphatFilmWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutRoughWallFunctionFvPatchScalarField_H
#define compressibleMutRoughWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
          Class alphatFilmWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatFilmWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- B Coefficient (default = 5.5)
        scalar B_;

        //- y+ value for laminar -> turbulent transition (default = 11.05)
        scalar yPlusCrit_;

        //- Turbulent Cmu coefficient (default = 0.09)
        scalar Cmu_;

        //- Von-Karman constant (default = 0.41)
        scalar kappa_;

        //- Turbulent Prandtl number (default = 0.85)
        scalar Prt_;


public:

    //- Runtime type information
    TypeName("alphatFilmWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatFilmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatFilmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatFilmWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatFilmWallFunctionFvPatchScalarField
        (
            const alphatFilmWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatFilmWallFunctionFvPatchScalarField
        (
            const alphatFilmWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatFilmWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphatFilmWallFunctionFvPatchScalarField
        (
            const alphatFilmWallFunctionFvPatchScalarField&,
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
                new alphatFilmWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

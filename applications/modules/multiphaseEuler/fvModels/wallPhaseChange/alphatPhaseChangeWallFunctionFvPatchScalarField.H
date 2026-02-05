/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::alphatPhaseChangeWallFunctionFvPatchScalarField

Description
    This boundary condition provides a thermal wall function or the turbulent
    thermal diffusivity for walls on which a wall phase change fvModel is in
    operation. The presence of this condition is used by the phase change
    fvModel to determine on which walls phase change takes place.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchType>
    {
        type            alphatPhaseChangeWallFunction;
        libs            ("libmultiphaseEulerFvModels.so");
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::fv::wallPhaseChange
    Foam::fv::wallBoiling
    Foam::fv::wallCondensation

SourceFiles
    alphatPhaseChangeWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatPhaseChangeWallFunctionFvPatchScalarField_H
#define alphatPhaseChangeWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "wallPhaseChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class alphatPhaseChangeWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatPhaseChangeWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
private:

    // Private Member Data

        //- List of the corresponding wall phase change models
        mutable UPtrList<const fv::wallPhaseChange> models_;


    // Private Member Functions

        //- Access the list of the corresponding wall phase change models
        const UPtrList<const fv::wallPhaseChange>& models() const;


public:

    //- Runtime type information
    TypeName("alphatPhaseChangeWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatPhaseChangeWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeWallFunctionFvPatchScalarField&,
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
                new alphatPhaseChangeWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

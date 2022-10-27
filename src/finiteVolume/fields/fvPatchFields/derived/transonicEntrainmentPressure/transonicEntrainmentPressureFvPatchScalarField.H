/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::transonicEntrainmentPressureFvPatchScalarField

Description
    This boundary condition provides an entrainment condition for pressure
    including support for supersonic jets exiting the domain.

Usage
    \table
        Property     | Description                 | Required | Default value
        rho          | Density field name          | no       | rho
        psi          | Compressibility field name  | no       | psi
        gamma        | Heat capacity ratio (cp/Cv) | yes      |
        Mb           | Blending Mach No            | no       | 0.5
        phi          | Flux field name             | no       | phi
        p0           | Reference pressure          | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            transonicEntrainmentPressure;
        gamma           1.4;
        p0              uniform 1e5;
    }
    \endverbatim

See also
    Foam::entrainmentPressureFvPatchScalarField
    Foam::mixedFvPatchField

SourceFiles
    transonicEntrainmentPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef transonicEntrainmentPressureFvPatchScalarField_H
#define transonicEntrainmentPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class transonicEntrainmentPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class transonicEntrainmentPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of density field field, default = "rho"
        const word rhoName_;

        //- Name of compressibility field, default = "psi"
        const word psiName_;

        //- Name of the flux field
        const word phiName_;

        //- Heat capacity ratio
        const scalar gamma_;

        //-
        const scalar Mb_;

        //- Reference pressure
        scalarField p0_;


public:

    //- Runtime type information
    TypeName("transonicEntrainmentPressure");


    // Constructors

        //- Construct from patch and internal field
        transonicEntrainmentPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transonicEntrainmentPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  transonicEntrainmentPressureFvPatchScalarField onto a new patch
        transonicEntrainmentPressureFvPatchScalarField
        (
            const transonicEntrainmentPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        transonicEntrainmentPressureFvPatchScalarField
        (
            const transonicEntrainmentPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        transonicEntrainmentPressureFvPatchScalarField
        (
            const transonicEntrainmentPressureFvPatchScalarField&,
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
                new transonicEntrainmentPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

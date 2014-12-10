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
    Foam::syringePressureFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a pressure condition, obtained from a
    zero-D model of the cylinder of a syringe.

    The syringe cylinder is defined by its initial volume, piston area and
    velocity profile specified by regions of constant acceleration, speed
    and deceleration.  The gas in the cylinder is described by its initial
    pressure and compressibility which is assumed constant, i.e. isothermal
    expansion/compression.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        Ap           | syringe piston area [m2] | yes        |
        Sp           | syringe piston speed [m/s] | yes      |
        VsI          | initial syringe volume [m3] | yes     |
        tas          | start of piston acceleration [s] | yes |
        tae          | end of piston acceleration [s] | yes  |
        tds          | start of piston deceleration [s] | yes |
        tde          | end of piston deceleration [s] | yes  |
        psI          | initial syringe pressure [Pa] | yes   |
        psi          | gas compressibility [m2/s2] | yes     |
        ams          | added (or removed) gas mass [kg] | yes |
    \endtable

    Example of the BC specification:
    \verbatim
    myPatch
    {
        type            syringePressure;
        Ap              1.388e-6;
        Sp              0.01;
        VsI             1.388e-8;
        tas             0.001;
        tae             0.002;
        tds             0.005;
        tde             0.006;
        psI             1e5;
        psi             1e-5;
        ams             0;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    syringePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef syringePressureFvPatchScalarField_H
#define syringePressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class syringePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class syringePressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- syringe piston area
        scalar Ap_;

        //- syringe piston speed
        scalar Sp_;

        //- Initial syringe volume
        scalar VsI_;

        //- Start of piston acceleration
        scalar tas_;

        //- End of piston acceleration
        scalar tae_;

        //- Start of piston deceleration
        scalar tds_;

        //- end of piston deceleration
        scalar tde_;

        //- Initial syringe pressure
        scalar psI_;

        //- Gas compressibility
        scalar psi_;

        //- Added gas mass
        scalar ams_;

        //- Added gas mass at previous time step
        scalar ams0_;

        //- Current time index used to store ms0_
        label curTimeIndex_;


    // Private Functions

        //- Return the volume of the syringe at time t
        scalar Vs(const scalar t) const;


public:

    //- Runtime type information
    TypeName("syringePressure");


    // Constructors

        //- Construct from patch and internal field
        syringePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        syringePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given syringePressureFvPatchScalarField
        //  onto a new patch
        syringePressureFvPatchScalarField
        (
            const syringePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new syringePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy
        syringePressureFvPatchScalarField
        (
            const syringePressureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        syringePressureFvPatchScalarField
        (
            const syringePressureFvPatchScalarField&,
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
                new syringePressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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

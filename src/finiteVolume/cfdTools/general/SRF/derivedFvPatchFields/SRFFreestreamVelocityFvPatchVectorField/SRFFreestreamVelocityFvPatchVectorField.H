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
    Foam::SRFFreestreamVelocityFvPatchVectorField

Description
    Freestream velocity condition to be used in conjunction with the single
    rotating frame (SRF) model (see: SRFModel class)

    Given the free stream velocity in the absolute frame, the condition
    applies the appropriate rotation transformation in time and space to
    determine the local velocity using:

        \f[
            U_p = cos(\theta)*U_{Inf} + sin(theta) (n^UInf) - U_{p,srf}
        \f]

    where
    \vartable
        U_p     = patch velocity [m/s]
        U_{Inf} = free stream velocity in the absolute frame [m/s]
        theta   = swept angle [rad]
        n       = axis direction of the SRF
        U_{p,srf} = SRF velocity of the patch
    \endvartable


Usage
    \table
        Property        | Description               | Required | Default value
        UInf            | freestream velocity       | yes      |
        relative        | UInf relative to the SRF? | no       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            SRFFreestreamVelocity;
        UInf            uniform (0 0 0);
        relative        no;
        value           uniform (0 0 0);    // initial value
    }
    \endverbatim

See also
    Foam::freestreamFvPatchField
    Foam::SRFVelocityFvPatchVectorField

SourceFiles
    SRFFreestreamVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef SRFFreestreamVelocityFvPatchVectorField_H
#define SRFFreestreamVelocityFvPatchVectorField_H

#include "inletOutletFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class SRFFreestreamVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class SRFFreestreamVelocityFvPatchVectorField
:
    public inletOutletFvPatchVectorField
{
    // Private Data

        //- Is the supplied inlet value relative to the SRF
        Switch relative_;

        //- Velocity of the free stream in the absolute frame [m/s]
        vector UInf_;


public:

    //- Runtime type information
    TypeName("SRFFreestreamVelocity");


    // Constructors

        //- Construct from patch and internal field
        SRFFreestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SRFFreestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given SRFFreestreamVelocityFvPatchVectorField
        //  onto a new patch
        SRFFreestreamVelocityFvPatchVectorField
        (
            const SRFFreestreamVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        SRFFreestreamVelocityFvPatchVectorField
        (
            const SRFFreestreamVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        SRFFreestreamVelocityFvPatchVectorField
        (
            const SRFFreestreamVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new SRFFreestreamVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the velocity at infinity
            const vector& UInf() const
            {
                return UInf_;
            }

            //- Return reference to the velocity at infinity to allow adjustment
            vector& UInf()
            {
                return UInf_;
            }


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

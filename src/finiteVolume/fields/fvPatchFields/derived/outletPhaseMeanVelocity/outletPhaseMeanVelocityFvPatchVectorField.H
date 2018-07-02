/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::outletPhaseMeanVelocityFvPatchVectorField

Description
    This boundary condition adjusts the velocity for the given phase to achieve
    the specified mean thus causing the phase-fraction to adjust according to
    the mass flow rate.

    Typical usage is as the outlet condition for a towing-tank ship simulation
    to maintain the outlet water level at the level as the inlet.

Usage
    \table
        Property     | Description             | Required    | Default value
        Umean        | mean velocity normal to the boundary [m/s] | yes |
        ramp         | ramping function for the mean flow speed | no | None
        alpha        | phase-fraction field    | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            outletPhaseMeanVelocity;
        Umean           1.2;
        alpha           alpha.water;
        value           uniform (1.2 0 0);
    }
    \endverbatim

See also
    Foam::mixedFvPatchField
    Foam::variableHeightFlowRateInletVelocityFvPatchVectorField

SourceFiles
    outletPhaseMeanVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef outletPhaseMeanVelocityFvPatchVectorField_H
#define outletPhaseMeanVelocityFvPatchVectorField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
    Class outletPhaseMeanVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class outletPhaseMeanVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        scalar Umean_;

        //- Ramp for the mean flow rate
        autoPtr<Function1<scalar>> ramp_;

        //- Name of the phase-fraction field
        word alphaName_;


public:

   //- Runtime type information
   TypeName("outletPhaseMeanVelocity");


   // Constructors

        //- Construct from patch and internal field
        outletPhaseMeanVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletPhaseMeanVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  outletPhaseMeanVelocityFvPatchVectorField
        //  onto a new patch
        outletPhaseMeanVelocityFvPatchVectorField
        (
            const outletPhaseMeanVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletPhaseMeanVelocityFvPatchVectorField
        (
            const outletPhaseMeanVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new outletPhaseMeanVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletPhaseMeanVelocityFvPatchVectorField
        (
            const outletPhaseMeanVelocityFvPatchVectorField&,
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
                new outletPhaseMeanVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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

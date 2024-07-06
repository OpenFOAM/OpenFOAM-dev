/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::swirlFlowRateInletVelocityFvPatchVectorField

Description
    Velocity inlet boundary condition creating a normal velocity field to match
    the specified mass or volumetric flow rate, and radial and tangential
    velocity fields specified by functions of time and radius or by a given
    angular speed.

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s

Usage
    \table
        Property     | Description             | Required    | Default value
        origin       | Origin of rotation      | no          | patch centre
        axis         | Axis of rotation        | no          | patch normal
        massFlowRate | Mass flow rate [kg/s]   | no          |
        volumetricFlowRate | Volumetric flow rate [m^3/s] | no |
        rho          | Density field name      | no          | rho
        rhoInlet     | Inlet density           | no          |
        radialVelocity | Radial velocity [m/s] | yes         |
        omega        | Angular velocity of the frame | no |
        tangentialVelocity | Tangential velocity [m/s] | no  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            swirlFlowRateInletVelocity;
        origin          (0 0 0);
        axis            (0 0 1);
        volumetricFlowRate constant 0.2;
        radialVelocity  constant 10;
        omega           constant 100 [rpm];
    }
    \endverbatim

    The \c volumetricFlowRate or \c massFlowRate and the \c omega entries are
    \c Function1 of time, see Foam::Function1s. The \c radialVelocity and
    \c tangentialVelocity entries are \c Function2 of time and radius,
    see Foam::Function2s.

    Note:
      - \c rhoInlet is required for the case of a mass flow rate, where the
        density field is not available at start-up
      - The value is positive into the domain (as an inlet)
      - May not work correctly for transonic inlets
      - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s
    Foam::Function2s

SourceFiles
    swirlFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef swirlFlowRateInletVelocityFvPatchVectorField_H
#define swirlFlowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function2.H"
#include "omega1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
        Class swirlFlowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class swirlFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Origin of the rotation
        const vector origin_;

        //- Axis of the rotation
        const vector axis_;

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Is the flow-rate volumetric?
        bool volumetric_;

        //- Name of the density field used to normalise the mass flux
        word rhoName_;

        //- Rho initialisation value (for start; if value not supplied)
        const scalar rhoInlet_;

        //- Radial velocity
        autoPtr<Function2<scalar>> radialVelocity_;

        //- Angular velocity of the frame
        autoPtr<Function1s::omega> omega_;

        //- Tangential velocity
        autoPtr<Function2<scalar>> tangentialVelocity_;


    // Private Member Functions

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues(const RhoType& rho);


public:

   //- Runtime type information
   TypeName("swirlFlowRateInletVelocity");


   // Constructors

        //- Construct from patch, internal field and dictionary
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&,
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
                new swirlFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::flowRateInletVelocityFvPatchVectorField

Description
    Velocity inlet boundary condition creating a velocity field with
    optionally specified profile normal to the patch adjusted to match the
    specified mass flow rate, volumetric flow rate or mean velocity.

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
        massFlowRate | Mass flow rate [kg/s]   | no          |
        volumetricFlowRate | Volumetric flow rate [m^3/s]| no |
        meanVelocity | Mean velocity [m/s]| no |
        profile      | Velocity profile        | no          |
        rho          | Density field name      | no          | rho
        rhoInlet     | Inlet density           | no          |
        alpha        | Volume fraction field name | no       |
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                flowRateInletVelocity;
        volumetricFlowRate  0.2;
        profile             laminarBL;
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
     \verbatim
    <patchName>
    {
        type                flowRateInletVelocity;
        massFlowRate        0.2;
        profile             turbulentBL;
        rho                 rho;
        rhoInlet            1.0;
    }
    \endverbatim

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                flowRateInletVelocity;
        meanVelocity        5;
        profile             turbulentBL;
    }
    \endverbatim

    The \c volumetricFlowRate, \c massFlowRate or \c meanVelocity entries are
    \c Function1 of time, see Foam::Function1s.

    The \c profile entry is a \c Function1 of the normalised distance to the
    wall.  Any suitable Foam::Function1s can be used including
    Foam::Function1s::codedFunction1 but Foam::Function1s::laminarBL and
    Foam::Function1s::turbulentBL have been created specifically for this
    purpose and are likely to be appropriate for most cases.

    Note:
      - \c rhoInlet is required for the case of a mass flow rate, where the
        density field is not available at start-up
      - The value is positive into the domain (as an inlet)
      - May not work correctly for transonic inlets
      - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s::laminarBL
    Foam::Function1s::turbulentBL
    Foam::Function1s
    Foam::flowRateOutletVelocityFvPatchVectorField

SourceFiles
    flowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletVelocityFvPatchVectorField_H
#define flowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class flowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class flowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Optional velocity profile
        autoPtr<Function1<scalar>> profile_;

        //- Is the flow-rate mean velocity?
        bool meanVelocity_;

        //- Is the flow-rate volumetric?
        bool volumetric_;

        //- Name of the density field used to normalise the mass flux
        word rhoName_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoInlet_;

        //- Name of the volume fraction field used to convert to interstitial
        //  velocity
        word alphaName_;

        //- Normalised distance to the wall
        scalarField y_;

        //- Total area of patch
        scalar area_;


    // Private Member Functions

        //- Set the wall distance field y_
        void setWallDist();

        //- Return the normalised velocity profile
        virtual tmp<scalarField> profile();

        //- Update the patch values given the appropriate volume fraction and
        //  density types and values
        template<class ScaleType, class AlphaType, class RhoType>
        void updateValues
        (
            const ScaleType& scale,
            const AlphaType& alpha,
            const RhoType& rho
        );

        //- Update the patch values given the appropriate volume fraction type
        //  and values
        template<class AlphaType>
        void updateValues(const AlphaType& alpha);

        //- Return whether the flow rate can be evaluated on this case (i.e.,
        //  is not a non-parallel processor case)
        bool canEvaluate();


public:

   //- Runtime type information
   TypeName("flowRateInletVelocity");


   // Constructors

        //- Construct from patch, internal field and dictionary
        flowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&,
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
                new flowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchVectorField&,
                const fvPatchFieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchVectorField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

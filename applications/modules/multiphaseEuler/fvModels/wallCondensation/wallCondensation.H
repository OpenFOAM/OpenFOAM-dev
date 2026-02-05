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
    Foam::fv::wallCondensation

Description
    Model for mass diffusion limited wall condensation between two phases on
    the surface of a number of wall patches. This model is based on an
    implementation described in Syrjanen, J., Hovi, V. (2023). The current
    implementation supports a single condensing specie.

    References:
    \verbatim
        Syrjanen, J., Hovi, V. (2023).
        Validation of a surface condensation model for PWR containment analysis.
        Research Report VTT-R-00884-23,
        VTT Technical Research Centre of Finland Ltd.
    \endverbatim

Usage
    Example usage:
    \verbatim
    wallCondensation
    {
        type            wallCondensation;
        libs            ("libmultiphaseEulerFvModels.so");

        // Note: Order is important. This model is one-way. It turns vapour
        // into liquid. The phases should be specified in this order.
        phases          (steam water);

        // The specie that is condensing
        specie  H2O;

        energySemiImplicit yes;
        specieSemiImplicit no;

        saturationPressure
        {
            type            ArdenBuck;
        }
    }
    \endverbatim

    In addition to the above fvModel specification, wall patches on which
    condensation is to be calculated should have an
    alphatPhaseChangeWallFunction boundary condition applied to the turbulent
    thermal diffusivity field.

See also
    Foam::alphatPhaseChangeWallFunctionFvPatchScalarField

SourceFiles
    wallCondensation.C

\*---------------------------------------------------------------------------*/

#ifndef wallCondensation_H
#define wallCondensation_H

#include "wallPhaseChange.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class saturationPressureModel;

class wallCondensationPhaseChangeRateFvPatchScalarField;

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class wallCondensation Declaration
\*---------------------------------------------------------------------------*/

class wallCondensation
:
    public wallPhaseChange
{
private:

    // Private classes

        //- Struct to hold various properties
        struct properties;


    // Private Data

        //- Reference to the liquid phase
        const phaseModel& liquid_;

        //- Reference to the vapour phase
        const phaseModel& vapour_;

        //- Reference to the liquid turbulent thermal diffusivity
        const volScalarField& alphatLiquid_;

        //- Reference to the vapour turbulent thermal diffusivity
        const volScalarField& alphatVapour_;

        //- Reference to the field associated with the pressure equation
        const volScalarField& p_rgh_;

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Turbulent Schmidt number
        scalar Sct_;

        //- The saturation curve
        autoPtr<saturationPressureModel> saturationModelPtr_;

        //- Counter for the evaluations of the pressure equation sources
        mutable label pressureEquationIndex_;

        //-  Whether or not to linearise the specie source
        bool specieSemiImplicit_;

        //- The phase change rate
        mutable volScalarField mDot_;

        //- The derivative of phase change rate
        mutable volScalarField mDotDy_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Correct the phase change rate
        void correctMDot() const;


public:

    //- Runtime type information
    TypeName("wallCondensation");


    // Constructors

        //- Construct from explicit source name and mesh
        wallCondensation
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Return a mask indicating whether phase change is occurring
            virtual const scalarField& active(const label) const;

            //- Access the turbulent thermal diffusivities for a patch
            virtual Pair<const scalarField&> alphats(const label) const;


        // Evaluation

            //- Return the fraction of the latent heat that is transferred into
            //  the second phase
            virtual tmp<DimensionedField<scalar, volMesh>> Lfraction() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Return the derivative of mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDotDy() const;

            //- Return the mass transfer rate for the given patch
            const wallCondensationPhaseChangeRateFvPatchScalarField& mDotPf
            (
                const label patchi
            ) const;

            //- Return the mass transfer rate for the given patch
            wallCondensationPhaseChangeRateFvPatchScalarField& mDotPfRef
            (
                const label patchi
            ) const;

            //- Use phaseChange's source functions
            using phaseChange::addSup;

            //- Override the pressure equation to update the mass transfer rate
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Override the species equations to linearise in the mass fraction
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& heOrYi,
                fvMatrix<scalar>& eqn
            ) const;


        //- Correct the fvModel
        virtual void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

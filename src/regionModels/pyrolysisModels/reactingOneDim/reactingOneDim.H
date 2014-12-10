/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::reactingOneDim

Description
    Reacting, 1-D pyrolysis model

SourceFiles
    reactingOneDim.C

\*---------------------------------------------------------------------------*/

#ifndef reactingOneDim_H
#define reactingOneDim_H

#include "pyrolysisModel.H"
#include "basicSolidChemistryModel.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{


/*---------------------------------------------------------------------------*\
                      Class reactingOneDim Declaration
\*---------------------------------------------------------------------------*/

class reactingOneDim
:
    public pyrolysisModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        reactingOneDim(const reactingOneDim&);

        //- Disallow default bitwise assignment
        void operator=(const reactingOneDim&);

        //- Read model controls
        void readReactingOneDimControls();


protected:

    // Protected data

        //- Reference to the solid chemistry model
        autoPtr<basicSolidChemistryModel> solidChemistry_;

        //- Reference to solid thermo
        solidReactionThermo& solidThermo_;

        //- Pointer to radiation model
        autoPtr<radiation::radiationModel> radiation_;


        // Reference to solid thermo properties

            //- Density [kg/m3]
            volScalarField rho_;

            //- List of solid components
            PtrList<volScalarField>& Ys_;

            // Non-const access to enthalpy
            volScalarField& h_;


        // Solution parameters

            //- Number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Maximum diffussivity
            scalar maxDiff_;

            //- Minimum delta for combustion
            scalar minimumDelta_;


        // Fields

            //- Total gas mass flux to the primary region [kg/m2/s]
            surfaceScalarField phiGas_;

            //- Sensible enthalpy gas flux [J/m2/s]
            volScalarField phiHsGas_;

            //- Heat release [J/s/m3]
            volScalarField chemistrySh_;


        // Source term fields

            //- Coupled region radiative heat flux [W/m2]
            //  Requires user to input mapping info for coupled patches
            //volScalarField QrCoupled_;

            //- In depth radiative heat flux [W/m2]
            volScalarField Qr_;


        // Checks

            //- Cumulative lost mass of the condensed phase [kg]
            dimensionedScalar lostSolidMass_;

            //- Cumulative mass generation of the gas phase [kg]
            dimensionedScalar addedGasMass_;

            //- Total mass gas flux at the pyrolysing walls [kg/s]
            scalar totalGasMassFlux_;

            //- Total heat release rate [J/s]
            dimensionedScalar totalHeatRR_;


        // Options

            //- Add gas enthalpy source term
            bool  gasHSource_;

            //- Add in depth radiation source term
            bool  QrHSource_;

            //- Use chemistry solvers (ode or sequential)
            bool useChemistrySolvers_;


    // Protected member functions

        //- Read control parameters from dictionary
        bool read();

        //- Read control parameters from dict
        bool read(const dictionary& dict);

        //- Update submodels
        void updateFields();

        //- Update/move mesh based on change in mass
        void updateMesh(const scalarField& mass0);

        //- Update radiative flux in pyrolysis region
        void updateQr();

        //- Update enthalpy flux for pyrolysis gases
        void updatePhiGas();

        //- Mass check
        void calculateMassTransfer();


        // Equations

            //- Solve continuity equation
            void solveContinuity();

            //- Solve energy
            void solveEnergy();

            //- Solve solid species mass conservation
            void solveSpeciesMass();


public:

    //- Runtime type information
    TypeName("reactingOneDim");


    // Constructors

        //- Construct from type name and mesh
        reactingOneDim
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name, mesh and dictionary
        reactingOneDim
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );




    //- Destructor
    virtual ~reactingOneDim();


    // Member Functions

        // Access

            //- Fields

                //- Return const density [Kg/m3]
                const volScalarField& rho() const;

                //- Return const temperature [K]
                virtual const volScalarField& T() const;

                //- Return specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const;

                //- Return the region absorptivity [1/m]
                virtual tmp<volScalarField> kappaRad() const;

                //- Return the region thermal conductivity [W/m/k]
                virtual tmp<volScalarField> kappa() const;

                //- Return the total gas mass flux to primary region [kg/m2/s]
                virtual const surfaceScalarField& phiGas() const;


        // Solution parameters

            //- Return the number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Return max diffusivity allowed in the solid
            virtual scalar maxDiff() const;


        // Helper functions

            //- External hook to add mass to the primary region
            virtual scalar addMassSources
            (
                const label patchI,            // patchI on primary region
                const label faceI              // faceI of patchI
            );

            //- Mean diffusion number of the solid region
            virtual scalar solidRegionDiffNo() const;


       // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the pyrolysis equations
            virtual void evolveRegion();


       // I-O

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactingOneDimI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

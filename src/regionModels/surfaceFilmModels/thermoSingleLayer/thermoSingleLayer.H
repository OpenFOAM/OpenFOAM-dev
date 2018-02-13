/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::thermoSingleLayer

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as Cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15K.  This is clearly non-conservative unless the heat-capacity
    is constant and should be rewritten to use the standard thermodynamics
    packages.

SourceFiles
    thermoSingleLayer.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayer_H
#define thermoSingleLayer_H

#include "kinematicSingleLayer.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// Forward declaration of classes
class filmViscosityModel;
class heatTransferModel;
class phaseChangeModel;
class filmRadiationModel;

/*---------------------------------------------------------------------------*\
                      Class thermoSingleLayer Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayer
:
    public kinematicSingleLayer
{
    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayer(const thermoSingleLayer&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayer&);

        //- Return boundary types for sensible enthalpy field
        wordList hsBoundaryTypes();


protected:

    // Protected data

        // Thermo properties

            //- Reference to the SLGThermo
            const SLGThermo& thermo_;


            // Fields

                //- Specific heat capacity [J/kg/K]
                volScalarField Cp_;

                //- Thermal conductivity [W/m/K]
                volScalarField kappa_;

                //- Temperature - mean [K]
                volScalarField T_;

                //- Temperature - surface [K]
                volScalarField Ts_;

                //- Temperature - wall [K]
                volScalarField Tw_;

                //- Sensible enthalpy [J/kg]
                volScalarField hs_;


            // Transfer fields - to the primary region

                //- Film energy transfer
                volScalarField primaryEnergyTrans_;


        //- Threshold film thickness beyond which the film is considered 'wet'
        scalar deltaWet_;


        // Hyprophilic/phobic properties

            //- Activation flag
            bool hydrophilic_;

            //- Length scale applied to deltaWet_ to determine when a wet
            //  surface becomes dry, typically 0.5
            scalar hydrophilicDryScale_;

            //- Length scale applied to deltaWet_ to determine when a dry
            //  surface becomes wet, typically 0.001
            scalar hydrophilicWetScale_;


        // Source term fields

            // Film region - registered to the film region mesh
            // Note: need boundary value mapped from primary region, and then
            // pushed into the patch internal field

                //- Energy [J/m2/s]
                volScalarField hsSp_;


            // Primary region - registered to the primary region mesh
            // Internal use only - not read-in

                //- Energy [J/m2/s]
                volScalarField hsSpPrimary_;


        // Fields mapped from primary region - registered to the film region
        // Note: need both boundary AND patch internal fields to be mapped

            //- Temperature [K]
            volScalarField TPrimary_;

            //- List of specie mass fractions [0-1]
            PtrList<volScalarField> YPrimary_;


        // Sub-models

            //- Viscosity model
            autoPtr<filmViscosityModel> viscosity_;

            //- Heat transfer coefficient between film surface and primary
            //  region [W/m2/K]
            autoPtr<heatTransferModel> htcs_;

            //- Heat transfer coefficient between wall and film [W/m2/K]
            autoPtr<heatTransferModel> htcw_;

            //- Phase change
            autoPtr<phaseChangeModel> phaseChange_;

            //- Radiation
            autoPtr<filmRadiationModel> radiation_;


        // Limits

            //- Minimum temperature limit (optional)
            scalar Tmin_;

            //- Maximum temperature limit (optional)
            scalar Tmax_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();

        //- Correct the thermo fields
        virtual void correctThermoFields();

        //- Correct sensible enthalpy for mapped temperature fields
        virtual void correctHsForMappedT();

        //- Correct the film surface and wall temperatures
        virtual void updateSurfaceTemperatures();

        //- Reset source term fields
        virtual void resetPrimaryRegionSourceTerms();

        //- Transfer thermo fields from the primary region to the film region
        virtual void transferPrimaryRegionThermoFields();

        //- Transfer source fields from the primary region to the film region
        virtual void transferPrimaryRegionSourceFields();

        //- Correct film coverage field
        virtual void correctAlpha();

        //- Update the film sub-models
        virtual void updateSubmodels();

        //- Return the wall/surface heat transfer term for the enthalpy equation
        virtual tmp<fvScalarMatrix> q(volScalarField& hs) const;


        // Equations

            //- Solve energy equation
            virtual void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermoSingleLayer");


    // Constructors

        //- Construct from components
        thermoSingleLayer
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType,
            const bool readFields = true
        );


    //- Destructor
    virtual ~thermoSingleLayer();


    // Member Functions

        // Thermo properties

            //- Return const reference to the SLGThermo object
            inline const SLGThermo& thermo() const;


            // Fields

                //- Return the film specific heat capacity [J/kg/K]
                virtual const volScalarField& Cp() const;

                //- Return the film thermal conductivity [W/m/K]
                virtual const volScalarField& kappa() const;

                //- Return the film mean temperature [K]
                virtual const volScalarField& T() const;

                //- Return the film surface temperature [K]
                virtual const volScalarField& Ts() const;

                //- Return the film wall temperature [K]
                virtual const volScalarField& Tw() const;

                //- Return the film sensible enthalpy [J/kg]
                virtual const volScalarField& hs() const;


            // Helper functions

                //- Return sensible enthalpy as a function of temperature
                //  for a patch
                inline tmp<scalarField> hs
                (
                    const scalarField& T,
                    const label patchi
                ) const;

                //- Return sensible enthalpy as a function of temperature
                inline tmp<volScalarField> hs
                (
                    const volScalarField& T
                ) const;

                //- Return temperature as a function of sensible enthalpy
                inline tmp<volScalarField> T
                (
                    const volScalarField& hs
                ) const;


         // Source fields (read/write access)

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchi,            // patchi on primary region
                const label facei,             // facei of patchi
                const scalar massSource,       // [kg]
                const vector& momentumSource,  // [kg.m/s] (tangential momentum)
                const scalar pressureSource,   // [kg.m/s] (normal momentum)
                const scalar energySource      // [J]
            );


        // Source term fields

            // Film region

                //- Energy [J/m2/s]
                inline const volScalarField& hsSp() const;


            // Primary region

                //- Energy [J/m2/s]
                inline const volScalarField& hsSpPrimary() const;


        // Fields mapped from the primary region

            //- Temperature [K]
            inline const volScalarField& TPrimary() const;

            //- Specie mass fractions [0-1]
            inline const PtrList<volScalarField>& YPrimary() const;



        // Sub-models

            //- Return const access to the (surface) heat transfer model
            inline const heatTransferModel& htcs() const;

            //- Return const access to the (wall) heat transfer model
            inline const heatTransferModel& htcw() const;

            //- Return const access to the phase change model
            inline const phaseChangeModel& phaseChange() const;

            //- Return const access to the radiation model
            inline const filmRadiationModel& radiation() const;


        // Derived fields (calculated on-the-fly)

            //- Return the convective heat energy from film to wall
            inline tmp<scalarField> Qconvw(const label patchi) const;

            //- Return the convective heat energy from primary region to film
            inline tmp<scalarField> Qconvp(const label patchi) const;


        // Evolution

            //- Pre-evolve film hook
            virtual void preEvolveRegion();

            //- Evolve the film equations
            virtual void evolveRegion();


        // Source fields

            // Mapped into primary region

                //- Return total mass source - Eulerian phase only
                virtual tmp<volScalarField::Internal> Srho() const;

                //- Return mass source for specie i - Eulerian phase only
                virtual tmp<volScalarField::Internal> Srho
                (
                    const label i
                ) const;

                //- Return enthalpy source - Eulerian phase only
                virtual tmp<volScalarField::Internal> Sh() const;


       // I-O

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoSingleLayerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

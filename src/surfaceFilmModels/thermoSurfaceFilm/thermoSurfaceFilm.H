/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::regionModels::thermoSurfaceFilm

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as Cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15K.  This is clearly non-conservative unless the heat-capacity
    is constant and should be rewritten to use the standard thermodynamics
    packages.

SourceFiles
    thermoSurfaceFilm.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSurfaceFilm_H
#define thermoSurfaceFilm_H

#include "momentumSurfaceFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fluidThermo;

namespace regionModels
{

// Forward declaration of classes
namespace surfaceFilmSubModels
{
    class viscosityModel;
    class heatTransferModel;
    class phaseChangeModel;
    class radiationModel;
}

/*---------------------------------------------------------------------------*\
                      Class thermoSurfaceFilm Declaration
\*---------------------------------------------------------------------------*/

class thermoSurfaceFilm
:
    public momentumSurfaceFilm
{

protected:

    // Protected data

        // Thermo properties

            //- Reference to the primary region thermo
            const fluidThermo& primaryThermo_;

            // Transfer fields - to the primary region

                //- Film energy transfer
                volScalarField primaryEnergyTrans_;


        //- Threshold film thickness beyond which the film is considered 'wet'
        scalar deltaWet_;


        // Hydrophilic/phobic properties

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
                volScalarField::Internal hSp_;


            // Primary region - registered to the primary region mesh
            // Internal use only - not read-in

                //- Energy [J/m2/s]
                volScalarField hSpPrimary_;


        // Fields mapped from primary region - registered to the film region
        // Note: need both boundary AND patch internal fields to be mapped

            //- Temperature [K]
            volScalarField TPrimary_;

            //- List of specie mass fractions [0-1]
            PtrList<volScalarField> YPrimary_;


        // Sub-models

            //- Heat transfer coefficient between film surface and primary
            //  region [W/m^2/K]
            autoPtr<surfaceFilmSubModels::heatTransferModel> htcs_;

            //- Heat transfer coefficient between wall and film [W/m^2/K]
            autoPtr<surfaceFilmSubModels::heatTransferModel> htcw_;

            //- Phase change
            autoPtr<surfaceFilmSubModels::phaseChangeModel> phaseChange_;

            //- Radiation
            autoPtr<surfaceFilmSubModels::radiationModel> radiation_;


        // Limits

            //- Minimum temperature limit (optional)
            scalar Tmin_;

            //- Maximum temperature limit (optional)
            scalar Tmax_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();

        //- Reset source term fields
        virtual void resetPrimaryRegionSourceTerms();

        //- Transfer thermo fields from the primary region to the film region
        virtual void transferPrimaryRegionThermoFields();

        //- Transfer source fields from the primary region to the film region
        virtual void transferPrimaryRegionSourceFields();

        //- Correct film coverage field
        virtual void correctCoverage();

        //- Update the film sub-models
        virtual void updateSubmodels();

        //- Return the wall/surface heat transfer term for the enthalpy equation
        virtual tmp<fvScalarMatrix> q(volScalarField& h) const;


        // Equations

            //- Solve energy equation
            virtual void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermoSurfaceFilm");


    // Constructors

        //- Construct from components
        thermoSurfaceFilm
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType,
            const bool readFields = true
        );

        //- Disallow default bitwise copy construction
        thermoSurfaceFilm(const thermoSurfaceFilm&) = delete;


    //- Destructor
    virtual ~thermoSurfaceFilm();


    // Member Functions

        // Thermo properties

            //- Return const reference to the primary region thermo object
            inline const fluidThermo& primaryThermo() const;


        // Derived Fields

            //- Return the film surface temperature [K]
            //  Currently this is assumed to be equal to
            //  the film mean temperature
            virtual tmp<volScalarField::Internal> Ts() const;

            //- Return the film wall temperature [K]
            virtual tmp<volScalarField::Internal> Tw() const;


         // Source fields (read/write access)

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchi,            // patchi on primary region
                const label facei,             // facei of patchi
                const scalar massSource,       // [kg]
                const vector& momentumSource,  // [kg m/s] (tangential momentum)
                const scalar pressureSource,   // [kg m/s] (normal momentum)
                const scalar energySource      // [J]
            );


        // Source term fields

            // Film region

                //- Energy [J/m2/s]
                inline const volScalarField::Internal& hSp() const;


            // Primary region

                //- Energy [J/m2/s]
                inline const volScalarField& hSpPrimary() const;


        // Fields mapped from the primary region

            //- Temperature [K]
            inline const volScalarField& TPrimary() const;

            //- Specie mass fractions [0-1]
            inline const PtrList<volScalarField>& YPrimary() const;



        // Sub-models

            //- Return const access to the (surface) heat transfer model
            inline const surfaceFilmSubModels::heatTransferModel& htcs() const;

            //- Return const access to the (wall) heat transfer model
            inline const surfaceFilmSubModels::heatTransferModel& htcw() const;

            //- Return const access to the phase change model
            inline const surfaceFilmSubModels::phaseChangeModel&
                phaseChange() const;

            //- Return const access to the radiation model
            inline const surfaceFilmSubModels::radiationModel&
                radiation() const;


        // Evolution

            //- Pre-evolve film hook
            virtual void preEvolveRegion();

            //- Evolve the film equations
            virtual void evolveRegion();


        // Source fields

            // Mapped into primary region

                //- Return mass source for specie i - Eulerian phase only
                virtual tmp<volScalarField::Internal> SYi
                (
                    const label i
                ) const;

                //- Return enthalpy source - Eulerian phase only
                virtual tmp<volScalarField::Internal> Sh() const;


       // I-O

            //- Provide some feedback
            virtual void info();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermoSurfaceFilm&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoSurfaceFilmI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

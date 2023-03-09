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
    Foam::thermoSurfaceFilm

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
namespace surfaceFilmModels
{
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

            //- Reference to the primary film thermo
            const fluidThermo& primaryThermo_;

            // Transfer fields - to the primary film

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

            // Film film - registered to the film film mesh
            // Note: need boundary value mapped from primary film, and then
            // pushed into the patch internal field

                //- Energy [J/m2/s]
                volScalarField::Internal hSp_;


            // Primary film - registered to the primary film mesh
            // Internal use only - not read-in

                //- Energy [J/m2/s]
                volScalarField hSpPrimary_;


        // Fields mapped from primary film - registered to the film film
        // Note: need both boundary AND patch internal fields to be mapped

            //- Temperature [K]
            volScalarField TPrimary_;

            //- List of specie mass fractions [0-1]
            PtrList<volScalarField> YPrimary_;


        // Sub-models

            //- Heat transfer coefficient between film surface and primary
            //  film [W/m^2/K]
            autoPtr<surfaceFilmModels::heatTransferModel> htcs_;

            //- Heat transfer coefficient between wall and film [W/m^2/K]
            autoPtr<surfaceFilmModels::heatTransferModel> htcw_;

            //- Phase change
            autoPtr<surfaceFilmModels::phaseChangeModel> phaseChange_;

            //- Radiation
            autoPtr<surfaceFilmModels::radiationModel> radiation_;


        // Limits

            //- Minimum temperature limit (optional)
            scalar Tmin_;

            //- Maximum temperature limit (optional)
            scalar Tmax_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();

        //- Reset source term fields
        virtual void resetPrimaryFilmSourceTerms();

        //- Transfer thermo fields from the primary film to the film film
        virtual void transferPrimaryFilmThermoFields();

        //- Transfer source fields from the primary film to the film film
        virtual void transferPrimaryFilmSourceFields();

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
            const fvMesh& primaryMesh,
            const dimensionedVector& g,
            const word& filmType,
            const bool readFields = true
        );

        //- Disallow default bitwise copy construction
        thermoSurfaceFilm(const thermoSurfaceFilm&) = delete;


    //- Destructor
    virtual ~thermoSurfaceFilm();


    // Member Functions

        // Thermo properties

            //- Return const reference to the primary film thermo object
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
                const label patchi,            // patchi on primary film
                const label facei,             // facei of patchi
                const scalar massSource,       // [kg]
                const vector& momentumSource,  // [kg m/s] (tangential momentum)
                const scalar pressureSource,   // [kg m/s] (normal momentum)
                const scalar energySource      // [J]
            );


        // Source term fields

            // Film film

                //- Energy [J/m2/s]
                inline const volScalarField::Internal& hSp() const;


            // Primary film

                //- Energy [J/m2/s]
                inline const volScalarField& hSpPrimary() const;


        // Fields mapped from the primary film

            //- Temperature [K]
            inline const volScalarField& TPrimary() const;

            //- Specie mass fractions [0-1]
            inline const PtrList<volScalarField>& YPrimary() const;



        // Sub-models

            //- Return const access to the (surface) heat transfer model
            inline const surfaceFilmModels::heatTransferModel& htcs() const;

            //- Return const access to the (wall) heat transfer model
            inline const surfaceFilmModels::heatTransferModel& htcw() const;

            //- Return const access to the phase change model
            inline const surfaceFilmModels::phaseChangeModel&
                phaseChange() const;

            //- Return const access to the radiation model
            inline const surfaceFilmModels::radiationModel&
                radiation() const;


        // Evolution

            //- Pre-evolve film hook
            virtual void preEvolveFilm();

            //- Evolve the film equations
            virtual void evolveFilm();


        // Source fields

            // Mapped into primary film

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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoSurfaceFilmI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

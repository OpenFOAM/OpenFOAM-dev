/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::MovingPhaseModel

Description
    Class which represents a moving fluid phase. Holds the velocity, fluxes and
    turbulence model and can generate the momentum equation. The interface is
    quite restrictive as it also has to support an equivalent stationary model,
    which does not store motion fields or a turbulence model.

    Possible future extensions include separating the turbulent functionality
    into another layer.

See also
    StationaryPhaseModel

SourceFiles
    MovingPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef MovingPhaseModel_H
#define MovingPhaseModel_H

#include "phaseModel.H"
#include "PhaseThermophysicalTransportModel.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Trait for converting the ThermoModel's thermo type to the thermo type needed
// for the thermophysical transport model type; i.e., from rho-type thermo to
// fluid-type thermo.

template<class ThermoModel>
struct MovingPhaseModelTransportThermoModel;

template<>
struct MovingPhaseModelTransportThermoModel<rhoThermo>
{
    typedef fluidThermo type;
};

template<>
struct MovingPhaseModelTransportThermoModel<rhoReactionThermo>
{
    typedef fluidReactionThermo type;
};

/*---------------------------------------------------------------------------*\
                      Class MovingPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class MovingPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected typedefs

        //- Thermo type for the thermophysical transport model
        typedef
            typename MovingPhaseModelTransportThermoModel
            <
                typename BasePhaseModel::thermoModel
            >::type
            transportThermoModel;


    // Protected data

        //- Velocity field
        volVectorField U_;

        //- Flux
        surfaceScalarField phi_;

        //- Volumetric flux
        surfaceScalarField alphaPhi_;

        //- Mass flux
        surfaceScalarField alphaRhoPhi_;

        //- Face velocity field
        tmp<surfaceVectorField> Uf_;

        //- Lagrangian acceleration field (needed for virtual-mass)
        mutable tmp<volVectorField> DUDt_;

        //- Lagrangian acceleration field on the faces (needed for virtual-mass)
        mutable tmp<surfaceScalarField> DUDtf_;

        //- Dilatation rate
        tmp<volScalarField> divU_;

        //- Turbulence model
        autoPtr<phaseCompressible::momentumTransportModel> turbulence_;

        //- Thermophysical transport model
        autoPtr
        <
            PhaseThermophysicalTransportModel
            <
                phaseCompressible::momentumTransportModel,
                transportThermoModel
            >
        > thermophysicalTransport_;

        //- Continuity error
        volScalarField continuityError_;

        //- Kinetic Energy
        mutable tmp<volScalarField> K_;


private:

    // Private static member functions

        //- Calculate and return the flux field
        tmp<surfaceScalarField> phi(const volVectorField& U) const;


public:

    // Constructors

        MovingPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~MovingPhaseModel();


    // Member Functions

        //- Correct the phase properties other than the thermo and turbulence
        virtual void correct();

        //- Correct the continuity error
        virtual void correctContinuityError(const volScalarField& source);

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport e.g. alphat
        virtual void correctEnergyTransport();

        //- Correct the face velocity for moving meshes
        virtual void correctUf();


        // Momentum

            //- Return whether the phase is stationary
            virtual bool stationary() const;

            //- Return the momentum equation
            virtual tmp<fvVectorMatrix> UEqn();

            //- Return the momentum equation for the face-based algorithm
            virtual tmp<fvVectorMatrix> UfEqn();

            //- Return the velocity
            virtual tmp<volVectorField> U() const;

            //- Access the velocity
            virtual volVectorField& URef();

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef();

            //- Return the face velocity
            //  Required for moving mesh cases
            virtual tmp<surfaceVectorField> Uf() const;

            //- Access the face velocity
            //  Required for moving mesh cases
            virtual surfaceVectorField& UfRef();

            //- Return the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhiRef();

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef();

            //- Return the substantive acceleration
            virtual tmp<volVectorField> DUDt() const;

            //- Return the substantive acceleration on the faces
            virtual tmp<surfaceScalarField> DUDtf() const;

            //- Return the continuity error
            virtual tmp<volScalarField> continuityError() const;

            //- Return the phase kinetic energy
            virtual tmp<volScalarField> K() const;


        // Compressibility (variable density)

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual tmp<volScalarField> divU() const;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual void divU(tmp<volScalarField> divU);


        // Momentum transport

            //- Return the turbulent kinetic energy
            virtual tmp<volScalarField> k() const;

            //- Return the phase-pressure'
            //  (derivative of phase-pressure w.r.t. phase-fraction)
            virtual tmp<volScalarField> pPrime() const;


        // Thermophysical transport

            //- Return the effective thermal conductivity on a patch
            virtual tmp<scalarField> kappaEff(const label patchi) const;

            //- Return the source term for the energy equation
            virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

            //- Return the source term for the given specie mass-fraction
            //  equation
            virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MovingPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

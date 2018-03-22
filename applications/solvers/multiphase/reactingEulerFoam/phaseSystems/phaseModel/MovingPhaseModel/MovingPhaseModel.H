/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    turbulence model. Provides access to the turbulent quantities.

    Possible future extensions include separating the turbulent fuctionality
    into another layer. It should also be possible to replace this layer with a
    stationary phase model, in order to model packed beds or simple porous
    media. This would probably require extra functionality, such as returning
    the inputs into the general pressure equation (A, HbyA, etc ...).

    Note that this class does not return the turbulence model, it just provides
    indirect access to the turbulent data. This is so a layer without
    turbulence modelling (such as a stationary model) could be substituted.

SourceFiles
    MovingPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef MovingPhaseModel_H
#define MovingPhaseModel_H

#include "phaseModel.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class MovingPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class MovingPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Velocity field
        volVectorField U_;

        //- Flux
        surfaceScalarField phi_;

        //- Volumetric flux
        surfaceScalarField alphaPhi_;

        //- Mass flux
        surfaceScalarField alphaRhoPhi_;

        //- Lagrangian acceleration field (needed for virtual-mass)
        mutable tmp<volVectorField> DUDt_;

        //- Lagrangian acceleration field on the faces (needed for virtual-mass)
        mutable tmp<surfaceScalarField> DUDtf_;

        //- Dilatation rate
        tmp<volScalarField> divU_;

        //- Turbulence model
        autoPtr<phaseCompressibleTurbulenceModel> turbulence_;

        //- Continuity error due to the flow
        volScalarField continuityErrorFlow_;

        //- Continuity error due to any sources
        volScalarField continuityErrorSources_;

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
            const label index
        );


    //- Destructor
    virtual ~MovingPhaseModel();


    // Member Functions

        //- Correct the phase properties other than the thermo and turbulence
        virtual void correct();

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport e.g. alphat
        virtual void correctEnergyTransport();


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

            //- Return the continuity error due to the flow field
            virtual tmp<volScalarField> continuityErrorFlow() const;

            //- Return the continuity error due to any sources
            virtual tmp<volScalarField> continuityErrorSources() const;

            //- Return the phase kinetic energy
            virtual tmp<volScalarField> K() const;


        // Compressibility (variable density)

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual tmp<volScalarField> divU() const;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual void divU(tmp<volScalarField> divU);


        // Transport (prevents compiler warnings)

            //- Return the effective thermal conductivity
            using BasePhaseModel::kappaEff;

            //- Return the effective thermal conductivity for enthalpy
            using BasePhaseModel::alphaEff;


        // Turbulence

            //- Return the turbulent dynamic viscosity
            virtual tmp<volScalarField> mut() const;

            //- Return the effective dynamic viscosity
            virtual tmp<volScalarField> muEff() const;

            //- Return the turbulent kinematic viscosity
            virtual tmp<volScalarField> nut() const;

            //- Return the effective kinematic viscosity
            virtual tmp<volScalarField> nuEff() const;

            //- Return the effective thermal conductivity
            virtual tmp<volScalarField> kappaEff() const;

            //- Return the effective thermal conductivity on a patch
            virtual tmp<scalarField> kappaEff(const label patchi) const;

            //- Return the effective thermal diffusivity for enthalpy
            virtual tmp<volScalarField> alphaEff() const;

            //- Return the effective thermal conductivity for enthalpy on a
            //  patch
            virtual tmp<scalarField> alphaEff(const label patchi) const;

            //- Return the turbulent kinetic energy
            virtual tmp<volScalarField> k() const;

            //- Return the phase-pressure'
            //  (derivative of phase-pressure w.r.t. phase-fraction)
            virtual tmp<volScalarField> pPrime() const;
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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class TransportModel>
class PhaseCompressibleTurbulenceModel;

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class MovingPhaseModel
:
    public BasePhaseModel
{
    // Private data

        //- Velocity field
        volVectorField U_;

        //- Flux
        surfaceScalarField phi_;

        //- Volumetric flux
        surfaceScalarField alphaPhi_;

        //- Mass flux
        surfaceScalarField alphaRhoPhi_;

        //- Lagrangian acceleration field (needed for virtual-mass)
        volVectorField DUDt_;

        //- Turbulence model
        autoPtr<PhaseCompressibleTurbulenceModel<phaseModel> > turbulence_;

        //- Continuity error
        volScalarField continuityError_;

        //- Phase diffusivity divided by the momentum coefficient.
        //  Used for implicit treatment of the phase pressure and dispersion
        tmp<surfaceScalarField> DbyA_;


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

        //- Return the momentum equation
        virtual tmp<fvVectorMatrix> UEqn();


        // Implicit phase pressure and dispersion support

            //- Return the phase diffusivity divided by the momentum coefficient
            virtual const surfaceScalarField& DbyA() const;

            //- Set the phase diffusivity divided by the momentum coefficient
            virtual void DbyA(const tmp<surfaceScalarField>& DbyA);


        // Momentum

            //- Constant access the velocity
            virtual tmp<volVectorField> U() const;

            //- Access the velocity
            virtual volVectorField& U();

            //- Return the substantive acceleration
            virtual tmp<volVectorField> DUDt() const;

            //- Constant access the continuity error
            virtual tmp<volScalarField> continuityError() const;

            //- Constant access the volumetric flux
            virtual tmp<surfaceScalarField> phi() const;

            //- Access the volumetric flux
            virtual surfaceScalarField& phi();

            //- Constant access the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhi();

            //- Constant access the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhi();


        // Turbulence

            //- Return the turbulence model
            virtual const PhaseCompressibleTurbulenceModel<phaseModel>&
                turbulence() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MovingPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

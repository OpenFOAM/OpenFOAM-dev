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
    Foam::phaseModel

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "transportModel.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class twoPhaseSystem;
class diameterModel;

template<class Phase>
class PhaseCompressibleTurbulenceModel;


/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField,
    public transportModel
{
    // Private data

        //- Reference to the twoPhaseSystem to which this phase belongs
        const twoPhaseSystem& fluid_;

        //- Name of phase
        word name_;

        dictionary phaseDict_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Thermophysical properties
        autoPtr<rhoThermo> thermo_;

        //- Velocity
        volVectorField U_;

        //- Volumetric flux of the phase
        surfaceScalarField alphaPhi_;

        //- Mass flux of the phase
        surfaceScalarField alphaRhoPhi_;

        //- Volumetric flux of the phase
        autoPtr<surfaceScalarField> phiPtr_;

        //- Diameter model
        autoPtr<diameterModel> dPtr_;

        //- turbulence model
        autoPtr<PhaseCompressibleTurbulenceModel<phaseModel> > turbulence_;


public:

    // Constructors

        phaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& phaseProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        const word& name() const
        {
            return name_;
        }

        //- Return the twoPhaseSystem to which this phase belongs
        const twoPhaseSystem& fluid() const
        {
            return fluid_;
        }

        //- Return the other phase in this two-phase system
        const phaseModel& otherPhase() const;

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return the Sauter-mean diameter
        tmp<volScalarField> d() const;

        //- Return the turbulence model
        const PhaseCompressibleTurbulenceModel<phaseModel>&
            turbulence() const;

        //- Return non-const access to the turbulence model
        //  for correction
        PhaseCompressibleTurbulenceModel<phaseModel>&
            turbulence();

        //- Return the thermophysical model
        const rhoThermo& thermo() const
        {
            return thermo_();
        }

        //- Return non-const access to the thermophysical model
        //  for correction
        rhoThermo& thermo()
        {
            return thermo_();
        }

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return thermo_->nu();
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return thermo_->nu(patchi);
        }

        //- Return the laminar dynamic viscosity
        tmp<volScalarField> mu() const
        {
            return thermo_->mu();
        }

        //- Return the laminar dynamic viscosity for patch
        tmp<scalarField> mu(const label patchi) const
        {
            return thermo_->mu(patchi);
        }

        //- Return the thermal conductivity
        tmp<volScalarField> kappa() const
        {
            return thermo_->kappa();
        }

        //- Return the specific heat capacity
        tmp<volScalarField> Cp() const
        {
            return thermo_->Cp();
        }

        //- Return the density
        const volScalarField& rho() const
        {
            return thermo_->rho();
        }

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return non-const access to the velocity
        //  Used in the momentum equation
        volVectorField& U()
        {
            return U_;
        }

        //- Return the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        //- Return non-const access to the volumetric flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- Return the volumetric flux of the phase
        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        //- Return non-const access to the volumetric flux of the phase
        surfaceScalarField& alphaPhi()
        {
            return alphaPhi_;
        }

        //- Return the mass flux of the phase
        const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Return non-const access to the mass flux of the phase
        surfaceScalarField& alphaRhoPhi()
        {
            return alphaRhoPhi_;
        }

        //- Correct the phase properties
        //  other than the thermodynamics and turbulence
        //  which have special treatment
        void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);

        //- Dummy Read for transportModel
        virtual bool read()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

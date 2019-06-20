/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
#include "fvMatricesFwd.H"
#include "rhoThermo.H"
#include "phaseCompressibleTurbulenceModelFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;
class diameterModel;

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
    // Private Data

        //- Reference to the phaseSystem to which this phase belongs
        const phaseSystem& fluid_;

        //- Name of phase
        word name_;

        //- Index of phase
        label index_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Diameter model
        autoPtr<diameterModel> diameterModel_;


public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            phaseSystem,
            (
                const phaseSystem& fluid,
                const word& phaseName,
                const label index
            ),
            (fluid, phaseName, index)
        );


    // Constructors

        phaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );

        //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const phaseSystem& fluid_;
            mutable label indexCounter_;

        public:

            iNew
            (
                const phaseSystem& fluid
            )
            :
                fluid_(fluid),
                indexCounter_(-1)
            {}

            autoPtr<phaseModel> operator()(Istream& is) const
            {
                indexCounter_++;
                return autoPtr<phaseModel>
                (
                    phaseModel::New(fluid_, word(is), indexCounter_)
                );
            }
        };


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        const word& name() const;

        //- Return the name of the phase for use as the keyword in PtrDictionary
        const word& keyword() const;

        //- Return the index of the phase
        label index() const;

        //- Return the system to which this phase belongs
        const phaseSystem& fluid() const;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const;

        //- Return the maximum phase-fraction (e.g. packing limit)
        scalar alphaMax() const;

        //- Return the Sauter-mean diameter
        tmp<volScalarField> d() const;

        //- Return const-reference to diameterModel of the phase
        const autoPtr<diameterModel>& dPtr() const;

        //- Correct the phase properties
        virtual void correct();

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport
        virtual void correctEnergyTransport();

        //- Ensure that the flux at inflow/outflow BCs is preserved
        void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

        //- Read phase properties dictionary
        virtual bool read();


        // Compressibility (variable density)

            //- Return true if the phase is compressible otherwise false
            virtual bool compressible() const = 0;

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual tmp<volScalarField> divU() const = 0;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual void divU(tmp<volScalarField> divU) = 0;


        // Thermo

            //- Return whether the phase is isothermal
            virtual bool isothermal() const = 0;

            //- Return the enthalpy equation
            virtual tmp<fvScalarMatrix> heEqn() = 0;

            //- Return the thermophysical model
            virtual const rhoThermo& thermo() const = 0;

            //- Access the thermophysical model
            virtual rhoThermo& thermoRef() = 0;

            //- Return the density field
            virtual tmp<volScalarField> rho() const = 0;


        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const = 0;

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi) = 0;

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const = 0;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef() = 0;

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& YActive() const = 0;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& YActiveRef() = 0;


        // Momentum

            //- Return whether the phase is stationary
            virtual bool stationary() const = 0;

            //- Return the momentum equation
            virtual tmp<fvVectorMatrix> UEqn() = 0;

            //- Return the momentum equation for the face-based algorithm
            virtual tmp<fvVectorMatrix> UfEqn() = 0;

            //- Return the velocity
            virtual tmp<volVectorField> U() const = 0;

            //- Access the velocity
            virtual volVectorField& URef() = 0;

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const = 0;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef() = 0;

            //- Return the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const = 0;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhiRef() = 0;

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const = 0;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef() = 0;

            //- Return the substantive acceleration
            virtual tmp<volVectorField> DUDt() const = 0;

            //- Return the substantive acceleration on the faces
            virtual tmp<surfaceScalarField> DUDtf() const = 0;

            //- Return the continuity error
            virtual tmp<volScalarField> continuityError() const = 0;

            //- Return the continuity error due to the flow field
            virtual tmp<volScalarField> continuityErrorFlow() const = 0;

            //- Return the continuity error due to any sources
            virtual tmp<volScalarField> continuityErrorSources() const = 0;

            //- Return the phase kinetic energy
            virtual tmp<volScalarField> K() const = 0;


        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const = 0;

            //- Return the laminar dynamic viscosity on a patch
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const = 0;

            //- Return the laminar kinematic viscosity on a patch
            virtual tmp<scalarField> nu(const label patchi) const = 0;

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual tmp<volScalarField> alpha() const = 0;

            //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alpha(const label patchi) const = 0;

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual tmp<volScalarField> kappa() const = 0;

            //- Thermal diffusivity for temperature of mixture
            //  for patch [W/m/K]
            virtual tmp<scalarField> kappa(const label patchi) const = 0;

            //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const = 0;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const = 0;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;


        // Turbulence

            //- Return the turbulent dynamic viscosity
            virtual tmp<volScalarField> mut() const = 0;

            //- Return the effective dynamic viscosity
            virtual tmp<volScalarField> muEff() const = 0;

            //- Return the turbulent kinematic viscosity
            virtual tmp<volScalarField> nut() const = 0;

            //- Return the effective kinematic viscosity
            virtual tmp<volScalarField> nuEff() const = 0;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff() const = 0;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff(const label patchi) const = 0;

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff() const = 0;

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff(const label patchi) const = 0;

            //- Return the turbulent kinetic energy
            virtual tmp<volScalarField> k() const = 0;

            //- Return the phase-pressure'
            //  (derivative of phase-pressure w.r.t. phase-fraction)
            virtual tmp<volScalarField> pPrime() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

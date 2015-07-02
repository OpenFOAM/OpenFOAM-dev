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
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;
class diameterModel;

template<class TransportModel>
class PhaseCompressibleTurbulenceModel;

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
    // Private data

        //- Reference to the phaseSystem to which this phase belongs
        const phaseSystem& fluid_;

        //- Name of phase
        word name_;

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
                const word& phaseName
            ),
            (fluid, phaseName)
        );


    // Constructors

        phaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Selectors

        static autoPtr<phaseModel> New
        (
            const phaseSystem& fluid,
            const word& phaseName
        );


    // Member Functions

        //- Return the name of this phase
        const word& name() const;

        //- Return the system to which this phase belongs
        const phaseSystem& fluid() const;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const;

        //- Return the maximum phase-fraction (e.g. packing limit)
        scalar alphaMax() const;

        //- Return the Sauter-mean diameter
        tmp<volScalarField> d() const;

        //- Correct the phase properties
        virtual void correct();

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Return the momentum equation
        virtual tmp<fvVectorMatrix> UEqn() = 0;

        //- Return the enthalpy equation
        virtual tmp<fvScalarMatrix> heEqn() = 0;

        //- Return the species fraction equation
        virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi) = 0;

        //- Read phase properties dictionary
        virtual bool read();


        // Compressibility (variable density)

            //- Return true if the phase is compressible otherwise false
            virtual bool compressible() const;

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual const volScalarField& divU() const;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual void divU(const volScalarField& divU);

            //- Return the phase kinetic energy
            virtual const volScalarField& K() const;


        // Implicit phase pressure and dispersion support

            //- Return the phase diffusivity divided by the momentum coefficient
            virtual const surfaceScalarField& DbyA() const;

            //- Set the phase diffusivity divided by the momentum coefficient
            virtual void DbyA(const tmp<surfaceScalarField>& DbyA);


        // Thermo

            //- Return const access to the thermophysical model
            virtual const rhoThermo& thermo() const = 0;

            //- Return non-const access to the thermophysical model
            //  for correction
            virtual rhoThermo& thermo() = 0;

            //- Return the density field
            virtual tmp<volScalarField> rho() const = 0;

            //- Constant access the species mass fractions
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& Y() = 0;


        // Momentum

            //- Constant access the velocity
            virtual tmp<volVectorField> U() const = 0;

            //- Access the velocity
            virtual volVectorField& U() = 0;

            //- Return the substantive acceleration
            virtual tmp<volVectorField> DUDt() const = 0;

            //- Constant access the continuity error
            virtual tmp<volScalarField> continuityError() const = 0;

            //- Constant access the volumetric flux
            virtual tmp<surfaceScalarField> phi() const = 0;

            //- Access the volumetric flux
            virtual surfaceScalarField& phi() = 0;

            //- Constant access the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const = 0;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhi() = 0;

            //- Constant access the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const = 0;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhi() = 0;


        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const = 0;

            //- Return the laminar dynamic viscosity on a patch
            virtual tmp<scalarField> mu(const label patchI) const = 0;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const = 0;

            //- Return the laminar kinematic viscosity on a patch
            virtual tmp<scalarField> nu(const label patchI) const = 0;

            //- Return the laminar thermal conductivity
            virtual tmp<volScalarField> kappa() const = 0;

            //- Return the laminar thermal conductivity on a patch
            virtual tmp<scalarField> kappa(const label patchI) const = 0;


        // Turbulence

            //- Return the turbulence model
            virtual const PhaseCompressibleTurbulenceModel<phaseModel>&
                turbulence() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

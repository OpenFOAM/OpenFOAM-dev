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
    Foam::ThermoPhaseModel

Description
    Class which represents a phase with a thermodynamic model. Provides access
    to the thermodynamic variables. Note that the thermo model itself is not
    returned as this class could be substituted in the hierarchy for one which
    mirrors the functionality, but does not include a thermo model; an
    incompressible phase model, for example.

SourceFiles
    ThermoPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoPhaseModel_H
#define ThermoPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class rhoThermo;

/*---------------------------------------------------------------------------*\
                      Class ThermoPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class ThermoType>
class ThermoPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Thermophysical model
        autoPtr<ThermoType> thermo_;


public:

    // Constructors

        ThermoPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~ThermoPhaseModel();


    // Member Functions

        // Thermo

            //- Return whether the phase is compressible
            virtual bool compressible() const;

            //- Return the thermophysical model
            virtual const rhoThermo& thermo() const;

            //- Access the thermophysical model
            virtual rhoThermo& thermoRef();

            //- Return the density field
            virtual tmp<volScalarField> rho() const;


        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const;

            //- Return the laminar dynamic viscosity on a patch
            virtual tmp<scalarField> mu(const label patchi) const;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const;

            //- Return the laminar kinematic viscosity on a patch
            virtual tmp<scalarField> nu(const label patchi) const;

            //- Return the laminar thermal conductivity
            virtual tmp<volScalarField> kappa() const;

            //- Return the laminar thermal conductivity on a patch
            virtual tmp<scalarField> kappa(const label patchi) const;

            //- Return the laminar thermal conductivity, given the turbulent
            //  thermal diffusivity
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const;

            //- Return the laminar thermal conductivity on a patch, given the
            //  turbulent thermal diffusivity
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Return the thermal diffusivity for enthalpy
            virtual tmp<volScalarField> alpha() const;

            //- Return the thermal diffusivity for enthalpy on a patch
            virtual tmp<scalarField> alpha(const label patchi) const;

            //- Return the effective thermal diffusivity for enthalpy, given the
            //  turbulent thermal diffusivity
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Return the effective thermal diffusivity for enthalpy on a
            //  patch, given the turbulent thermal diffusivity
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;


        // Turbulence (prevents compiler warnings)

            //- Return the effective thermal conductivity
            using BasePhaseModel::kappaEff;

            //- Return the effective thermal conductivity for enthalpy
            using BasePhaseModel::alphaEff;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

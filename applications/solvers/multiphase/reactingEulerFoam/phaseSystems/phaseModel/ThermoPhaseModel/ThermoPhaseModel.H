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
    Foam::ThermoPhaseModel

Description
    Class which represents a phase with a thermodynamic model. Provides access
    to the thermodynamic variables. Note that the thermo model itself is not
    returned as this class could be substituted in the hierarcy for one which
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
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class ThermoType>
class ThermoPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Thermophysical model
        autoPtr<ThermoType> thermoPtr_;

        //- Thermophysical model
        ThermoType* thermo_;


public:

    // Constructors

        ThermoPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index,
            const bool createThermo = true
        );


    //- Destructor
    virtual ~ThermoPhaseModel();


    // Member Functions

        // Thermo

            //- Return const access to the thermophysical model
            virtual const rhoThermo& thermo() const;

            //- Return non-const access to the thermophysical model
            //  for correction
            virtual rhoThermo& thermo();

            //- Return the density field
            virtual tmp<volScalarField> rho() const;


        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const;

            //- Access the laminar dynamic viscosity
            virtual tmp<scalarField> mu(const label patchi) const;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const;

            //- Access the laminar kinematic viscosity
            virtual tmp<scalarField> nu(const label patchi) const;

            //- Return the laminar thermal conductivity
            virtual tmp<volScalarField> kappa() const;

            //- Access the laminar thermal conductivity
            virtual tmp<scalarField> kappa(const label patchi) const;

            //- Return the laminar thermal conductivity
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const;

            //- Access the laminar thermal conductivity
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Return the thermal diffusivity for enthalpy
            virtual tmp<volScalarField> alpha() const;

            //- Return the thermal diffusivity for enthalpy on a patch
            virtual tmp<scalarField> alpha(const label patchi) const;

            //- Return the thermal diffusivity for enthalpy
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Return the thermal diffusivity for enthalpy on a patch
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ThermoPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

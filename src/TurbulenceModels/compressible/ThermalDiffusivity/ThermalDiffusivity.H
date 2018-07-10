/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::ThermalDiffusivity

Description
    Templated wrapper class to provide compressible turbulence models
    thermal diffusivity based thermal transport.

SourceFiles
    ThermalDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef ThermalDiffusivity_H
#define ThermalDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class ThermalDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class ThermalDiffusivity
:
    public BasicTurbulenceModel
{

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef volScalarField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct
        ThermalDiffusivity
        (
            const word& type,
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& trasport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<ThermalDiffusivity> New
        (
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& trasportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


        //- Return a reference to the selected turbulence model
        static autoPtr<ThermalDiffusivity> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const transportModel& trasportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~ThermalDiffusivity()
    {}


    // Member Functions

        //- Return the thermal diffusivity for temperature [J/m/s/K]
        virtual tmp<volScalarField> kappa() const
        {
            return this->transport_.kappa();
        }

        //- Return the laminar thermal diffusivity for temperature on patch
        //  [J/m/s/K]
        virtual tmp<scalarField> kappa(const label patchi) const
        {
            return this->transport_.kappa(patchi);
        }

        //- Return the laminar thermal diffusivity for enthalpy [kg/m/s]
        virtual tmp<volScalarField> alpha() const
        {
            return this->transport_.alpha();
        }

        //- Return the laminar thermal diffusivity for enthalpy on patch
        //  [kg/m/s]
        virtual tmp<scalarField> alpha(const label patchi) const
        {
            return this->transport_.alpha(patchi);
        }

        //- Return the turbulent thermal diffusivity for enthalpy [kg/m/s]
        virtual tmp<volScalarField> alphat() const;

        //- Return the turbulent thermal diffusivity for enthalpy for a patch
        //  [kg/m/s]
        virtual tmp<scalarField> alphat(const label patchi) const;

        //- Return the effective turbulent thermal diffusivity for temperature
        //  [J/m/s/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return kappa();
        }

        //- Return the effective turbulent thermal diffusivity for temperature
        //  [J/m/s/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return kappa(patchi);
        }

        //- Return the effective turbulent thermal diffusivity for enthalpy
        //  [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            return alpha();
        }

        //- Return the effective turbulent thermal diffusivity for enthalpy
        //  for a patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return alpha(patchi);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermalDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

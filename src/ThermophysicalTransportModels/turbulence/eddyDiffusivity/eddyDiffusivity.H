/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::turbulenceThermophysicalTransportModels::eddyDiffusivity

Description
    Eddy-diffusivity based temperature gradient heat flux model
    for single specie RAS or LES of turbulent flow.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

Usage
    \verbatim
    LES
    {
        model           eddyDiffusivity;
        Prt             0.85;
    }
    \endverbatim

SourceFiles
    eddyDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDiffusivity_H
#define eddyDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                       Class eddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceThermophysicalTransportModel>
class eddyDiffusivity
:
    public TurbulenceThermophysicalTransportModel
{

protected:

    // Protected data

        // Model coefficients

            //- Turbulent Prandtl number []
            dimensionedScalar Prt_;

        // Fields

            //- Turbulent thermal diffusivity of enthalpy [kg/m/s]
            volScalarField alphat_;


    // Protected Member Functions

        virtual void correctAlphat();


public:

    typedef typename TurbulenceThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        TurbulenceThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename TurbulenceThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("eddyDiffusivity");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        eddyDiffusivity
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~eddyDiffusivity()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Turbulent thermal diffusivity for enthalpy [kg/m/s]
        virtual tmp<volScalarField> alphat() const
        {
            return alphat_;
        }

        //- Turbulent thermal diffusivity for enthalpy for a patch [kg/m/s]
        virtual tmp<scalarField> alphat(const label patchi) const
        {
            return alphat()().boundaryField()[patchi];
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return this->thermo().kappa() + this->thermo().Cp()*alphat();
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return
                this->thermo().kappa().boundaryField()[patchi]
              + this->thermo().Cp().boundaryField()[patchi]*alphat(patchi);
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            return this->thermo().alphahe() + alphat();
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return this->thermo().alphahe(patchi) + alphat(patchi);
        }

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Correct the eddyDiffusivity viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "eddyDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

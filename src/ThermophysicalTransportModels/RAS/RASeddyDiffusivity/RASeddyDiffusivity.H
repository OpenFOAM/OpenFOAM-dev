/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::RASThermophysicalTransportModels::eddyDiffusivity

Description
    Eddy-diffusivity based gradient heat flux model for RAS flow.

SourceFiles
    RASeddyDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef RASeddyDiffusivity_H
#define RASeddyDiffusivity_H

#include "RASThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASThermophysicalTransportModels
{

/*---------------------------------------------------------------------------* \
                           Class eddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class eddyDiffusivity
:
    public RASThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Prt_;

        // Fields

            volScalarField alphat_;


    // Protected Member Functions

        virtual void correctAlphat();


public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;


    //- Runtime type information
    TypeName("eddyDiffusivity");


    // Constructors

        //- Construct from components
        eddyDiffusivity(const momentumTransportModel& momentumTransport);


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
            return this->thermo().kappaEff(alphat());
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return this->thermo().kappaEff
            (
                alphat(patchi),
                patchi
            );
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            return this->thermo().alphaEff(alphat());
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return this->thermo().alphaEff
            (
                alphat(patchi),
                patchi
            );
        }

        //- Return the heat flux
        virtual tmp<volVectorField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Correct the eddyDiffusivity viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RASeddyDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

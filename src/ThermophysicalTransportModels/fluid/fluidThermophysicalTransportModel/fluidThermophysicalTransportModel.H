/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::fluidThermophysicalTransportModel

Description
    Abstract base class for fluid thermophysical transport models
    RAS, LES and laminar.

SourceFiles
    fluidThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermophysicalTransportModel_H
#define fluidThermophysicalTransportModel_H

#include "thermophysicalTransportModel.H"
#include "compressibleMomentumTransportModel.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fluidThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

class fluidThermophysicalTransportModel
:
    public thermophysicalTransportModel
{
protected:

    // Protected data

        const compressibleMomentumTransportModel& momentumTransportModel_;


public:

    // Constructors

        //- Construct from compressibleMomentumTransportModel
        fluidThermophysicalTransportModel
        (
            const compressibleMomentumTransportModel& momentumTransport
        );


    //- Destructor
    virtual ~fluidThermophysicalTransportModel()
    {}


    // Member Functions

        //- Access function to momentum transport model
        const compressibleMomentumTransportModel& momentumTransport() const
        {
            return momentumTransportModel_;
        }

        //- Access function to fluid thermophysical properties
        virtual const fluidThermo& thermo() const = 0;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const = 0;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const = 0;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const = 0;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const = 0;

        //- Return the patch heat flux correction [W/m^2]
        //  For isotropic or patch-aligned thermal conductivity qCorr is null
        virtual tmp<scalarField> qCorr(const label patchi) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const = 0;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

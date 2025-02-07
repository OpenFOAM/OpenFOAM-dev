/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Foam::ThermophysicalTransportPhaseModel

Description
    Class which makes thermophysical transport modelling available to derived
    classes

SourceFiles
    ThermophysicalTransportPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef ThermophysicalTransportPhaseModel_H
#define ThermophysicalTransportPhaseModel_H

#include "phaseModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "PhaseThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidThermo;
class rhoFluidThermo;
class fluidMulticomponentThermo;
class rhoFluidMulticomponentThermo;

// Trait for converting the ThermoModel's thermo type to the thermo type needed
// for the thermophysical transport model type; i.e., from rho-type thermo to
// fluid-type thermo.

template<class ThermoModel>
struct MovingPhaseModelTransportThermoModel;

template<>
struct MovingPhaseModelTransportThermoModel<rhoFluidThermo>
{
    typedef fluidThermo type;
};

template<>
struct MovingPhaseModelTransportThermoModel<rhoFluidMulticomponentThermo>
{
    typedef fluidMulticomponentThermo type;
};


/*---------------------------------------------------------------------------*\
                   Class ThermophysicalTransportPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class ThermophysicalTransportPhaseModel
:
    public BasePhaseModel
{
    // Private Member Data

        //- Thermophysical transport model type
        typedef
            PhaseThermophysicalTransportModel
            <
                phaseCompressible::momentumTransportModel,
                typename MovingPhaseModelTransportThermoModel
                <
                    typename BasePhaseModel::thermoModel
                >::type
            >
            thermophysicalTransportModel;

        //- Thermophysical transport model
        autoPtr<thermophysicalTransportModel> thermophysicalTransport_;


public:

    // Constructors

        ThermophysicalTransportPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~ThermophysicalTransportPhaseModel();


    // Member Functions

        //- Predict the energy transport e.g. alphat
        virtual void predictThermophysicalTransport();

        //- Correct the energy transport e.g. alphat
        virtual void correctThermophysicalTransport();

        //- Return the effective thermal conductivity on a patch
        virtual tmp<scalarField> kappaEff(const label patchi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the source term for a specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermophysicalTransportPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

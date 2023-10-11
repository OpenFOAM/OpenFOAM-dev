/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
    Foam::AnisothermalPhaseModel

Description
    Class which represents a phase for which the temperature (strictly energy)
    varies. Returns the energy equation and corrects the thermodynamic model.

SourceFiles
    AnisothermalPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef AnisothermalPhaseModel_H
#define AnisothermalPhaseModel_H

#include "phaseModel.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "PhaseThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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
                   Class AnisothermalPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class AnisothermalPhaseModel
:
    public BasePhaseModel
{
    // Private Member Data

        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;

        //- Thermo type for the thermophysical transport model
        typedef
            typename MovingPhaseModelTransportThermoModel
            <
                typename BasePhaseModel::thermoModel
            >::type
            transportThermoModel;

        //- Thermophysical transport model
        autoPtr
        <
            PhaseThermophysicalTransportModel
            <
                phaseCompressible::momentumTransportModel,
                transportThermoModel
            >
        > thermophysicalTransport_;


    // Private Member Functions

        //- Optionally filter the pressure work term as the phase-fraction -> 0
        tmp<volScalarField> filterPressureWork
        (
            const tmp<volScalarField>& pressureWork
        ) const;


public:

    // Constructors

        AnisothermalPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~AnisothermalPhaseModel();


    // Member Functions

        //- Return whether the phase is isothermal
        virtual bool isothermal() const;

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Predict the energy transport e.g. alphat
        virtual void predictThermophysicalTransport();

        //- Correct the energy transport e.g. alphat
        virtual void correctThermophysicalTransport();

        //- Return the effective thermal conductivity on a patch
        virtual tmp<scalarField> kappaEff(const label patchi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the source term for the given specie mass-fraction
        //  equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Return the enthalpy equation
        virtual tmp<fvScalarMatrix> heEqn();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AnisothermalPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

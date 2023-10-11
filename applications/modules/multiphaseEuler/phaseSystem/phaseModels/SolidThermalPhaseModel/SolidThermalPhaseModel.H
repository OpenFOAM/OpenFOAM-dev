/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::SolidThermalPhaseModel

Description
    Class which represents a solid stationary phase for which the temperature
    (strictly energy) varies. Returns the energy equation and corrects the
    thermodynamic model.

SourceFiles
    SolidThermalPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef SolidThermalPhaseModel_H
#define SolidThermalPhaseModel_H

#include "phaseModel.H"
#include "phaseSolidThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SolidThermalPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class SolidThermalPhaseModel
:
    public BasePhaseModel
{
    // Thermophysical transport

        //- Pointer to the solid thermophysical transport model
        autoPtr<phaseSolidThermophysicalTransportModel>
            thermophysicalTransport_;


public:

    // Constructors

        SolidThermalPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~SolidThermalPhaseModel();


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

        //- Return the enthalpy equation
        virtual tmp<fvScalarMatrix> heEqn();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolidThermalPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

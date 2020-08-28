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
    Foam::PhaseThermophysicalTransportModel

Description
    Templated base class for multiphase thermophysical transport models.

SourceFiles
    PhaseThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseThermophysicalTransportModel_H
#define PhaseThermophysicalTransportModel_H

#include "ThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class PhaseThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class MomentumTransportModel, class ThermoModel>
class PhaseThermophysicalTransportModel
:
    public ThermophysicalTransportModel<MomentumTransportModel, ThermoModel>
{

public:

    typedef MomentumTransportModel momentumTransportModel;
    typedef ThermoModel thermoModel;
    typedef typename momentumTransportModel::alphaField alphaField;


    // Constructors

        //- Construct
        PhaseThermophysicalTransportModel
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<PhaseThermophysicalTransportModel> New
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~PhaseThermophysicalTransportModel()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhaseThermophysicalTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

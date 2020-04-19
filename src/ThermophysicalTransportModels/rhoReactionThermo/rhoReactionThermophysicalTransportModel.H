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

Typedef
    Foam::compressible::thermophysicalTransportModel

Typedef
    Foam::compressible::laminarThermophysicalTransportModel

Typedef
    Foam::compressible::RASThermophysicalTransportModel

Typedef
    Foam::compressible::LESThermophysicalTransportModel

Description
    Typedefs for laminar, RAS and LES thermophysical transport models
    for rhoReactionThermo flow.

\*---------------------------------------------------------------------------*/

#ifndef rhoReactionThermophysicalTransportModel_H
#define rhoReactionThermophysicalTransportModel_H

#include "rhoReactionThermo.H"
#include "laminarThermophysicalTransportModel.H"
#include "RASThermophysicalTransportModel.H"
#include "LESThermophysicalTransportModel.H"
#include "fluidThermoMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ThermophysicalTransportModel
    <
        compressible::momentumTransportModel,
        rhoReactionThermo
    > rhoReactionThermophysicalTransportModel;

    typedef laminarThermophysicalTransportModel
    <
        rhoReactionThermophysicalTransportModel
    > laminarRhoReactionThermophysicalTransportModel;

    typedef RASThermophysicalTransportModel
    <
        rhoReactionThermophysicalTransportModel
    > RASRhoReactionThermophysicalTransportModel;

    typedef LESThermophysicalTransportModel
    <
        rhoReactionThermophysicalTransportModel
    > LESRhoReactionThermophysicalTransportModel;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

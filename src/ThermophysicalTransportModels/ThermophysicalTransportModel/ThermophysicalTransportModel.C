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

\*---------------------------------------------------------------------------*/

#include "ThermophysicalTransportModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class ThermoModel>
Foam::ThermophysicalTransportModel<MomentumTransportModel, ThermoModel>::
ThermophysicalTransportModel
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    thermophysicalTransportModel(momentumTransport),
    momentumTransport_(momentumTransport),
    thermo_(thermo)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class ThermoModel>
Foam::autoPtr
<
    Foam::ThermophysicalTransportModel<MomentumTransportModel, ThermoModel>
>
Foam::ThermophysicalTransportModel<MomentumTransportModel, ThermoModel>::New
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
{
    const word modelType
    (
        momentumTransport.lookup("simulationType")
    );

    Info<< "Selecting thermophysical transport type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown thermophysical transport type "
            << modelType << nl << nl
            << "Available types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ThermophysicalTransportModel>
    (
        cstrIter()(momentumTransport, thermo)
    );
}


// ************************************************************************* //

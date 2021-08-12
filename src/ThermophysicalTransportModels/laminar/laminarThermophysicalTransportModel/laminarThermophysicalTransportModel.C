/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "laminarThermophysicalTransportModel.H"
#include "unityLewisFourier.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
void Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::printCoeffs
(
    const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::laminarThermophysicalTransportModel
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    BasicThermophysicalTransportModel(momentumTransport, thermo),
    laminarDict_(this->subOrEmptyDict("laminar")),
    printCoeffs_(laminarDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(laminarDict_.optionalSubDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Foam::autoPtr
<
    Foam::laminarThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >
>
Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::New
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
{
    typeIOobject<IOdictionary> header
    (
        IOobject::groupName
        (
            thermophysicalTransportModel::typeName,
            momentumTransport.alphaRhoPhi().group()
        ),
        momentumTransport.time().constant(),
        momentumTransport.mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (header.headerOk())
    {
        IOdictionary modelDict(header);

        const word modelType(modelDict.subDict("laminar").lookup( "model"));

        Info<< "Selecting laminar thermophysical transport model "
            << modelType << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown laminar thermophysical transport model "
                << modelType << nl << nl
                << "Available models:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<laminarThermophysicalTransportModel>
        (
            cstrIter()(momentumTransport, thermo)
        );
    }
    else
    {
        Info<< "Selecting default laminar thermophysical transport model "
            << laminarThermophysicalTransportModels::unityLewisFourier<
               BasicThermophysicalTransportModel>::typeName << endl;

        return autoPtr<laminarThermophysicalTransportModel>
        (
            new laminarThermophysicalTransportModels::unityLewisFourier
            <
                BasicThermophysicalTransportModel
            >(momentumTransport, thermo)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::read()
{
    if (BasicThermophysicalTransportModel::read())
    {
        laminarDict_ <<= this->subDict("laminar");

        coeffDict_ <<= laminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicThermophysicalTransportModel>
void Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::correct()
{
    BasicThermophysicalTransportModel::correct();
}


// ************************************************************************* //

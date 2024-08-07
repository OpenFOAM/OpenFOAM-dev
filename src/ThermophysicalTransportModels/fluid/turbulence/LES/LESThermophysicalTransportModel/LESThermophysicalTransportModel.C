/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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

#include "LESThermophysicalTransportModel.H"
#include "unityLewisEddyDiffusivity.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Foam::LESThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::LESThermophysicalTransportModel
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    BasicThermophysicalTransportModel(momentumTransport, thermo)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Foam::autoPtr
<
    Foam::LESThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >
>
Foam::LESThermophysicalTransportModel
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

        const word modelType(modelDict.subDict("LES").lookup( "model"));

        Info<< "Selecting LES thermophysical transport model "
            << modelType << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown LES thermophysical transport model "
                << modelType << nl << nl
                << "Available models:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<LESThermophysicalTransportModel>
        (
            cstrIter()(momentumTransport, thermo)
        );
    }
    else
    {
        typedef
            turbulenceThermophysicalTransportModels::unityLewisEddyDiffusivity
            <
                LESThermophysicalTransportModel
                <
                    BasicThermophysicalTransportModel
                >
            > LESunityLewisEddyDiffusivity;

        Info<< "Selecting default LES thermophysical transport model "
            <<  LESunityLewisEddyDiffusivity::typeName << endl;

        return autoPtr<LESThermophysicalTransportModel>
        (
            new LESunityLewisEddyDiffusivity
            (
                LESunityLewisEddyDiffusivity::typeName,
                momentumTransport,
                thermo,
                true
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
const Foam::dictionary& Foam::LESThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::coeffDict() const
{
    return this->subOrEmptyDict("LES").optionalSubDict(type() + "Coeffs");
}


template<class BasicThermophysicalTransportModel>
bool Foam::LESThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::read()
{
    return BasicThermophysicalTransportModel::read();
}


template<class BasicThermophysicalTransportModel>
void Foam::LESThermophysicalTransportModel<BasicThermophysicalTransportModel>::
predict()
{
    BasicThermophysicalTransportModel::predict();
}


template<class BasicThermophysicalTransportModel>
void Foam::LESThermophysicalTransportModel<BasicThermophysicalTransportModel>::
correct()
{
    BasicThermophysicalTransportModel::correct();
}


// ************************************************************************* //

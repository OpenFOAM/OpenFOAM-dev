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

#include "RASThermophysicalTransportModel.H"
#include "unityLewisEddyDiffusivity.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
void Foam::RASThermophysicalTransportModel
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
Foam::RASThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::RASThermophysicalTransportModel
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    BasicThermophysicalTransportModel(momentumTransport, thermo),
    RASDict_(this->subOrEmptyDict("RAS")),
    printCoeffs_(RASDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(RASDict_.optionalSubDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Foam::autoPtr
<
    Foam::RASThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >
>
Foam::RASThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::New
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
{
    IOobject header
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

    if (header.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary modelDict(header);

        const word modelType(modelDict.subDict("RAS").lookup( "model"));

        Info<< "Selecting RAS thermophysical transport model "
            << modelType << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown RAS thermophysical transport model "
                << modelType << nl << nl
                << "Available models:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<RASThermophysicalTransportModel>
        (
            cstrIter()(momentumTransport, thermo)
        );
    }
    else
    {
        typedef
            turbulenceThermophysicalTransportModels::unityLewisEddyDiffusivity
            <
                RASThermophysicalTransportModel
                <
                    BasicThermophysicalTransportModel
                >
            > RASunityLewisEddyDiffusivity;

        Info<< "Selecting default RAS thermophysical transport model "
            <<  RASunityLewisEddyDiffusivity::typeName << endl;

        return autoPtr<RASThermophysicalTransportModel>
        (
            new RASunityLewisEddyDiffusivity
            (
                RASunityLewisEddyDiffusivity::typeName,
                momentumTransport,
                thermo,
                true
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool Foam::RASThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::read()
{
    if (BasicThermophysicalTransportModel::read())
    {
        RASDict_ <<= this->subDict("RAS");

        coeffDict_ <<= RASDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicThermophysicalTransportModel>
void Foam::RASThermophysicalTransportModel<BasicThermophysicalTransportModel>::
correct()
{
    BasicThermophysicalTransportModel::correct();
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

#include "turbulenceFields.H"
#include "incompressibleMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulenceFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        turbulenceFields,
        dictionary
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFields::compressibleField,
    7
>::names[] =
{
    "k",
    "epsilon",
    "omega",
    "nut",
    "nuEff",
    "kappaEff",
    "R"
};

const Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFields::compressibleField,
    7
> Foam::functionObjects::turbulenceFields::compressibleFieldNames_;

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFields::incompressibleField,
    6
>::names[] =
{
    "k",
    "epsilon",
    "omega",
    "nut",
    "nuEff",
    "R"
};

const Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFields::incompressibleField,
    6
> Foam::functionObjects::turbulenceFields::incompressibleFieldNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFields::turbulenceFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFields::~turbulenceFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceFields::read(const dictionary& dict)
{
    if (dict.found("field"))
    {
        fieldSet_.insert(word(dict.lookup("field")));
    }
    else
    {
        fieldSet_.insert(wordList(dict.lookup("fields")));
    }

    if (dict.lookupOrDefault<Switch>("prefix", false))
    {
        prefix_ = momentumTransportModel::typeName + ':';
    }
    else
    {
        prefix_ = word::null;
    }

    Info<< type() << " " << name() << ": ";
    if (fieldSet_.size())
    {
        Info<< "storing fields:" << nl;
        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            Info<< "    "
                << IOobject::groupName(prefix_ + iter.key(), phaseName_) << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::turbulenceFields::execute()
{
    if (obr_.foundType<fluidThermophysicalTransportModel>(phaseName_))
    {
        const fluidThermophysicalTransportModel& ttm =
            obr_.lookupType<fluidThermophysicalTransportModel>(phaseName_);

        const compressibleMomentumTransportModel& model =
            ttm.momentumTransport();

        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            const word& f = iter.key();
            switch (compressibleFieldNames_[f])
            {
                case compressibleField::k:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case compressibleField::epsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case compressibleField::omega:
                {
                    processField<scalar>(f, model.omega());
                    break;
                }
                case compressibleField::nut:
                {
                    processField<scalar>(f, model.nut());
                    break;
                }
                case compressibleField::nuEff:
                {
                    processField<scalar>(f, model.nuEff());
                    break;
                }
                case compressibleField::kappaEff:
                {
                    processField<scalar>(f, ttm.kappaEff());
                    break;
                }
                case compressibleField::R:
                {
                    processField<symmTensor>(f, model.sigma());
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << exit(FatalError);
                }
            }
        }
    }
    else if (obr_.foundType<compressibleMomentumTransportModel>(phaseName_))
    {
        const compressibleMomentumTransportModel& model =
            obr_.lookupType<compressibleMomentumTransportModel>(phaseName_);

        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            const word& f = iter.key();
            switch (compressibleFieldNames_[f])
            {
                case compressibleField::k:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case compressibleField::epsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case compressibleField::omega:
                {
                    processField<scalar>(f, model.omega());
                    break;
                }
                case compressibleField::nut:
                {
                    processField<scalar>(f, model.nut());
                    break;
                }
                case compressibleField::nuEff:
                {
                    processField<scalar>(f, model.nuEff());
                    break;
                }
                case compressibleField::R:
                {
                    processField<symmTensor>(f, model.sigma());
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << exit(FatalError);
                }
            }
        }
    }
    else if
    (
        obr_.foundType<incompressible::momentumTransportModel>(phaseName_)
    )
    {
        const incompressible::momentumTransportModel& model =
            obr_.lookupType<incompressible::momentumTransportModel>(phaseName_);

        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            const word& f = iter.key();
            switch (incompressibleFieldNames_[f])
            {
                case incompressibleField::k:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case incompressibleField::epsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case incompressibleField::omega:
                {
                    processField<scalar>(f, model.omega());
                    break;
                }
                case incompressibleField::nut:
                {
                    processField<scalar>(f, model.nut());
                    break;
                }
                case incompressibleField::nuEff:
                {
                    processField<scalar>(f, model.nuEff());
                    break;
                }
                case incompressibleField::R:
                {
                    processField<symmTensor>(f, model.sigma());
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << exit(FatalError);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Turbulence model not found in database, deactivating"
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::turbulenceFields::write()
{
    forAllConstIter(wordHashSet, fieldSet_, iter)
    {
        const word fieldName
        (
            IOobject::groupName(prefix_ + iter.key(), phaseName_)
        );
        writeObject(fieldName);
    }

    return true;
}


// ************************************************************************* //

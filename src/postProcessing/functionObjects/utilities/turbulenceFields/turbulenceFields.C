/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
#include "dictionary.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulenceFields, 0);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFields::compressibleField,
    8
>::names[] =
{
    "k",
    "epsilon",
    "mut",
    "muEff",
    "alphat",
    "alphaEff",
    "R",
    "devRhoReff"
};

const Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFields::compressibleField,
    8
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
    "nut",
    "nuEff",
    "R",
    "devReff"
};

const Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFields::incompressibleField,
    6
> Foam::functionObjects::turbulenceFields::incompressibleFieldNames_;

const Foam::word Foam::functionObjects::turbulenceFields::modelName
(
    Foam::turbulenceModel::propertiesName
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceFields::compressible()
{
    if (obr_.foundObject<compressible::turbulenceModel>(modelName))
    {
        return true;
    }
    else if (obr_.foundObject<incompressible::turbulenceModel>(modelName))
    {
        return false;
    }
    else
    {
        FatalErrorInFunction
            << "Turbulence model not found in database, deactivating"
            << exit(FatalError);
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFields::turbulenceFields
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    fieldSet_()
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    if
    (
       !obr.foundObject<compressible::turbulenceModel>(modelName)
    && !obr.foundObject<incompressible::turbulenceModel>(modelName)
    )
    {
        FatalErrorInFunction
            << "Cannot find turbulenceModel in objectRegistry"
            << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFields::~turbulenceFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::turbulenceFields::read(const dictionary& dict)
{
    fieldSet_.insert(wordList(dict.lookup("fields")));

    Info<< type() << " " << name_ << ": ";
    if (fieldSet_.size())
    {
        Info<< "storing fields:" << nl;
        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            Info<< "    " << modelName << ':' << iter.key() << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }
}


void Foam::functionObjects::turbulenceFields::execute()
{
    bool comp = compressible();

    if (comp)
    {
        const compressible::turbulenceModel& model =
            obr_.lookupObject<compressible::turbulenceModel>(modelName);

        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            const word& f = iter.key();
            switch (compressibleFieldNames_[f])
            {
                case cfK:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case cfEpsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case cfMut:
                {
                    processField<scalar>(f, model.mut());
                    break;
                }
                case cfMuEff:
                {
                    processField<scalar>(f, model.muEff());
                    break;
                }
                case cfAlphat:
                {
                    processField<scalar>(f, model.alphat());
                    break;
                }
                case cfAlphaEff:
                {
                    processField<scalar>(f, model.alphaEff());
                    break;
                }
                case cfR:
                {
                    processField<symmTensor>(f, model.R());
                    break;
                }
                case cfDevRhoReff:
                {
                    processField<symmTensor>(f, model.devRhoReff());
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }
    else
    {
        const incompressible::turbulenceModel& model =
            obr_.lookupObject<incompressible::turbulenceModel>(modelName);

        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            const word& f = iter.key();
            switch (incompressibleFieldNames_[f])
            {
                case ifK:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case ifEpsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case ifNut:
                {
                    processField<scalar>(f, model.nut());
                    break;
                }
                case ifNuEff:
                {
                    processField<scalar>(f, model.nuEff());
                    break;
                }
                case ifR:
                {
                    processField<symmTensor>(f, model.R());
                    break;
                }
                case ifDevReff:
                {
                    processField<symmTensor>(f, model.devReff());
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }
}


void Foam::functionObjects::turbulenceFields::end()
{
    execute();
}


void Foam::functionObjects::turbulenceFields::timeSet()
{}


void Foam::functionObjects::turbulenceFields::write()
{}


// ************************************************************************* //

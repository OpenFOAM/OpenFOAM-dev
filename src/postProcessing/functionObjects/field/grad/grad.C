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

#include "grad.H"
#include "volFields.H"
#include "dictionary.H"
#include "grad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(grad, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::grad::grad
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    fieldName_("undefined-fieldName"),
    resultName_("undefined-resultName")
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::grad::~grad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::grad::read(const dictionary& dict)
{
    dict.lookup("fieldName") >> fieldName_;
    dict.lookup("resultName") >> resultName_;

    if (resultName_ == "none")
    {
        resultName_ = "fvc::grad(" + fieldName_ + ")";
    }
}


void Foam::functionObjects::grad::execute()
{
    bool processed = false;

    calcGrad<scalar>(fieldName_, resultName_, processed);
    calcGrad<vector>(fieldName_, resultName_, processed);

    if (!processed)
    {
        WarningInFunction
            << "Unprocessed field " << fieldName_ << endl;
    }
}


void Foam::functionObjects::grad::end()
{
    execute();
}


void Foam::functionObjects::grad::timeSet()
{}


void Foam::functionObjects::grad::write()
{
    if (obr_.foundObject<regIOobject>(resultName_))
    {
        const regIOobject& field =
            obr_.lookupObject<regIOobject>(resultName_);

        Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << field.name() << nl << endl;

        field.write();
    }
}


// ************************************************************************* //

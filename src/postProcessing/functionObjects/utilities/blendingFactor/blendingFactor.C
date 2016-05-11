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

#include "blendingFactor.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blendingFactor, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blendingFactor::blendingFactor
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    phiName_("unknown-phiName"),
    fieldName_("unknown-fieldName")
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blendingFactor::~blendingFactor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::blendingFactor::read(const dictionary& dict)
{
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    dict.lookup("fieldName") >> fieldName_;
}


void Foam::functionObjects::blendingFactor::execute()
{
    calc<scalar>();
    calc<vector>();
}


void Foam::functionObjects::blendingFactor::end()
{
    execute();
}

void Foam::functionObjects::blendingFactor::timeSet()
{}


void Foam::functionObjects::blendingFactor::write()
{
    const word fieldName = "blendingFactor:" + fieldName_;

    const volScalarField& blendingFactor =
        obr_.lookupObject<volScalarField>(fieldName);

    Info<< type() << " " << name_ << " output:" << nl
        << "    writing field " << blendingFactor.name() << nl
        << endl;

    blendingFactor.write();
}


// ************************************************************************* //

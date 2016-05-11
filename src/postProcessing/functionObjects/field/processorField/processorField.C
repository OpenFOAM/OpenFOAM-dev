/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "processorField.H"
#include "dictionary.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(processorField, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::processorField::processorField
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr)
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volScalarField* procFieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                "processorID",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimless, 0.0)
        )
    );

    mesh.objectRegistry::store(procFieldPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::processorField::~processorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::processorField::read(const dictionary& dict)
{}


void Foam::functionObjects::processorField::execute()
{
    const volScalarField& procField =
        obr_.lookupObject<volScalarField>("processorID");

    const_cast<volScalarField&>(procField) ==
        dimensionedScalar("proci", dimless, Pstream::myProcNo());
}


void Foam::functionObjects::processorField::end()
{
    execute();
}


void Foam::functionObjects::processorField::timeSet()
{}


void Foam::functionObjects::processorField::write()
{
    const volScalarField& procField =
        obr_.lookupObject<volScalarField>("processorID");

    procField.write();
}


// ************************************************************************* //

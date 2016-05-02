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

#include "calcFvcDiv.H"
#include "volFields.H"
#include "dictionary.H"
#include "calcFvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(calcFvcDiv, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::functionObjects::calcFvcDiv::divField
(
    const word& divName,
    const dimensionSet& dims
)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    if (!mesh.foundObject<volScalarField>(divName))
    {
        volScalarField* divFieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    divName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dims/dimLength, 0.0)
            )
        );

        mesh.objectRegistry::store(divFieldPtr);
    }

    const volScalarField& field = mesh.lookupObject<volScalarField>(divName);

    return const_cast<volScalarField&>(field);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::calcFvcDiv::calcFvcDiv
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
    read(dict);
}


bool Foam::functionObjects::calcFvcDiv::viable
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
{
    // Construction is viable if the available mesh is an fvMesh
    return isA<fvMesh>(obr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::calcFvcDiv::~calcFvcDiv()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::calcFvcDiv::read(const dictionary& dict)
{
    dict.lookup("fieldName") >> fieldName_;
    dict.lookup("resultName") >> resultName_;

    if (resultName_ == "none")
    {
        resultName_ = "fvc::div(" + fieldName_ + ")";
    }
}


void Foam::functionObjects::calcFvcDiv::execute()
{
    bool processed = false;

    calcDiv<surfaceScalarField>(fieldName_, resultName_, processed);
    calcDiv<volVectorField>(fieldName_, resultName_, processed);

    if (!processed)
    {
        WarningInFunction
            << "Unprocessed field " << fieldName_ << endl;
    }
}


void Foam::functionObjects::calcFvcDiv::end()
{
    execute();
}


void Foam::functionObjects::calcFvcDiv::timeSet()
{}


void Foam::functionObjects::calcFvcDiv::write()
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

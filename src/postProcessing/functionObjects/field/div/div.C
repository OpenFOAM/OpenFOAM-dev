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

#include "div.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(div, 0);
    addToRunTimeSelectionTable(functionObject, div, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::functionObjects::div::divField
(
    const word& divName,
    const dimensionSet& dims
)
{
    if (!mesh_.foundObject<volScalarField>(divName))
    {
        volScalarField* divFieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    divName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dims/dimLength, 0.0)
            )
        );

        mesh_.objectRegistry::store(divFieldPtr);
    }

    const volScalarField& field = mesh_.lookupObject<volScalarField>(divName);

    return const_cast<volScalarField&>(field);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::div::div
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::div::~div()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::div::read(const dictionary& dict)
{
    dict.lookup("fieldName") >> fieldName_;
    dict.lookup("resultName") >> resultName_;

    if (resultName_ == "none")
    {
        resultName_ = "fvc::div(" + fieldName_ + ")";
    }

    return true;
}


bool Foam::functionObjects::div::execute(const bool postProcess)
{
    bool processed = false;

    calcDiv<surfaceScalarField>(fieldName_, resultName_, processed);
    calcDiv<volVectorField>(fieldName_, resultName_, processed);

    if (!processed)
    {
        WarningInFunction
            << "Unprocessed field " << fieldName_ << endl;
    }

    return true;
}


bool Foam::functionObjects::div::write(const bool postProcess)
{
    if (obr_.foundObject<regIOobject>(resultName_))
    {
        const regIOobject& field =
            obr_.lookupObject<regIOobject>(resultName_);

        Info<< type() << " " << name() << " output:" << nl
            << "    writing field " << field.name() << nl << endl;

        field.write();
    }

    return true;
}


// ************************************************************************* //

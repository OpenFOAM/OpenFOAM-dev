/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "turbulenceIntensity.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulenceIntensity, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        turbulenceIntensity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::turbulenceIntensity::writeFileHeader(const label i)
{
    writeHeader(file(), "I ()");
    writeCommented(file(), "Time");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceIntensity::turbulenceIntensity
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log)
{
    volScalarField* turbulenceIntensityPtr
    (
        new volScalarField
        (
            IOobject
            (
                "I",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0)
        )
    );

    mesh_.objectRegistry::store(turbulenceIntensityPtr);

    read(dict);
    resetName(typeName);
    resetLocalObjectName("I");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceIntensity::~turbulenceIntensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceIntensity::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    return true;
}


bool Foam::functionObjects::turbulenceIntensity::execute()
{
    volScalarField& turbulenceIntensity =
        mesh_.lookupObjectRef<volScalarField>("I");

    if (mesh_.foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        const turbulenceModel& turbModel = mesh_.lookupObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

        volScalarField uPrime(sqrt((2.0/3.0)*turbModel.k()));
        turbulenceIntensity =
            uPrime
           /max
            (
                max(uPrime, mag(turbModel.U())),
                dimensionedScalar("small", dimVelocity, small)
            );
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::turbulenceIntensity::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& turbulenceIntensity =
        mesh_.lookupObject<volScalarField>("I");

    const scalar minTurbulenceIntensity = min(turbulenceIntensity).value();
    const scalar maxTurbulenceIntensity = max(turbulenceIntensity).value();
    const scalar avgTurbulenceIntensity = turbulenceIntensity.average().value();

    if (Pstream::master())
    {
        Log << "    I : min = " << minTurbulenceIntensity
            << ", max = " << maxTurbulenceIntensity
            << ", average = " << avgTurbulenceIntensity << nl;

        writeTime(file());
        file()
            << tab << minTurbulenceIntensity
            << tab << maxTurbulenceIntensity
            << tab << avgTurbulenceIntensity
            << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //

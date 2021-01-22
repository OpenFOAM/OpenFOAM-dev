/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
#include "momentumTransportModel.H"
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
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceIntensity::~turbulenceIntensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceIntensity::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    resetName(typeName);
    resetLocalObjectName("I");

    return true;
}


bool Foam::functionObjects::turbulenceIntensity::execute()
{
    if
    (
        mesh_.foundObject<momentumTransportModel>
        (
            momentumTransportModel::typeName
        )
    )
    {
        const momentumTransportModel& turbModel =
            mesh_.lookupObject<momentumTransportModel>
            (
                momentumTransportModel::typeName
            );

        volScalarField uPrime(sqrt((2.0/3.0)*turbModel.k()));

        word name("I");

        return
            store
            (
                name,
                uPrime
               /max
                (
                    max(uPrime, mag(turbModel.U())),
                    dimensionedScalar(dimVelocity, small)
                )
            );
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);

        return false;
    }
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

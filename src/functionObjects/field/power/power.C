/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "power.H"
#include "momentumTransportModel.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(power, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        power,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * * //

const Foam::wordList Foam::functionObjects::power::fields_
{
    "stressUSf", "divStressU", "divStressDotU", "stressDDotGradU",
    "pUSf", "divPU", "gradPU", "pDivU",
    "tauUSf", "divTauU", "divTauDotU", "tauDDotGradU"
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::power::writeFileHeader(const label i)
{
    file().setf(ios_base::fixed, ios_base::floatfield);

    writeHeader
    (
        file(),
        "Domain-integrated power [x " + Foam::name(1.0/factor_) + "]"
    );
    writeCommented(file(), "Time");

    forAll(fields_, fieldi)
    {
        // Do not write surface fields (0, 4, 8)
        if ((fieldi % 4) != 0)
        {
            writeTabbed(file(), fields_[fieldi]);
        }
    }

    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::power::power
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    factor_(1e9)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::power::~power()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::power::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    dict.readIfPresent("factor", factor_);

    resetName(typeName);
    resetLocalObjectNames(fields_);

    return true;
}


bool Foam::functionObjects::power::execute()
{
    const momentumTransportModel& transport =
        mesh_.lookupType<momentumTransportModel>();

    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volVectorField& U = transport.U();
    const surfaceVectorField Uf(fvc::interpolate(U));

    // Deviatoric stress tensor (the spherical part is added to the pressure)
    const surfaceVectorField tauf(-transport.devTau()*mesh().magSf());

    // Pressure stress
    const surfaceVectorField pIf(fvc::interpolate(p)*mesh().Sf());

    // Full deviatoric stress tensor = tau - pI
    const surfaceVectorField stressf(tauf - pIf);

    store
    (
        fields_[3],
        store(fields_[1], fvc::div(store(fields_[0], stressf & Uf)))
      - store(fields_[2], fvc::div(stressf) & U)
    );

    store
    (
        fields_[7],
        store(fields_[5], fvc::div(store(fields_[4], -pIf & Uf)))
      - store(fields_[6], -(fvc::grad(p) & U))
    );

    store
    (
        fields_[11],
        store(fields_[9], fvc::div(store(fields_[8], tauf & Uf)))
      - store(fields_[10], fvc::div(tauf) & U)
    );

    return true;
}


bool Foam::functionObjects::power::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    forAll(fields_, fieldi)
    {
        if (mesh_.foundObject<volScalarField>(fields_[fieldi]))
        {
            const volScalarField& power =
                mesh_.lookupObject<volScalarField>(fields_[fieldi]);

            const scalar integratedPower =
                factor_*fvc::domainIntegrate(power).value();

            if (Pstream::master())
            {
                Log << "    " << fields_[fieldi] << ": "
                    << integratedPower << nl;

                file()<< tab << integratedPower;
            }
        }
    }

    if (Pstream::master())
    {
        Log << endl;
        file() << endl;
    }

    return true;
}


// ************************************************************************* //

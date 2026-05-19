/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "bXiQdot.H"
#include "ubRhoThermo.H"
#include "uRhoMulticomponentThermo.H"
#include "bRhoMulticomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bXiQdot, 0);
    addToRunTimeSelectionTable(functionObject, bXiQdot, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::bXiQdot::bXiQdot
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_)
{
    read(dict);
    resetLocalObjectName(typeName);

    // Request caching of the b-equation source
    obr_.cacheTemporary("bSource");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::bXiQdot::~bXiQdot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::bXiQdot::read(const dictionary& dict)
{
    functionObject::read(dict);
    writeLocalObjects::read(dict);

    return true;
}


bool Foam::functionObjects::bXiQdot::execute()
{
    const word uThermoName
    (
        IOobject::groupName
        (
            physicalProperties::typeName,
            ubRhoThermo::unburntPhaseName
        )
    );

    const word bThermoName
    (
        IOobject::groupName
        (
            physicalProperties::typeName,
            ubRhoThermo::burntPhaseName
        )
    );

    if
    (
        mesh_.foundObject<uRhoMulticomponentThermo>(uThermoName)
     && mesh_.foundObject<bRhoMulticomponentThermo>(bThermoName)
    )
    {
        const uRhoMulticomponentThermo& uThermo =
            mesh_.lookupObject<uRhoMulticomponentThermo>(uThermoName);

        const bRhoMulticomponentThermo& bThermo =
            mesh_.lookupObject<bRhoMulticomponentThermo>(bThermoName);

        const volScalarField::Internal& bSource
        (
            mesh_.lookupObject<volScalarField::Internal>("bSource")
        );

        store
        (
            typeName,
            bSource*(bThermo.hf()()() - uThermo.hf()()())
        );

        return true;
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find "
               "uRhoMulticomponentThermo or bRhoMulticomponentThermo "
               "in the database"
            << exit(FatalError);

        return false;
    }
}


bool Foam::functionObjects::bXiQdot::write()
{
    return writeLocalObjects::write();
}


// ************************************************************************* //

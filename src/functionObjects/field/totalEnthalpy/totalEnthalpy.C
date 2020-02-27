/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "totalEnthalpy.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(totalEnthalpy, 0);
    addToRunTimeSelectionTable(functionObject, totalEnthalpy, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::totalEnthalpy::totalEnthalpy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_, false),
    phaseName_(word::null)
{
    read(dict);
    resetLocalObjectName(IOobject::groupName("Ha", phaseName_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::totalEnthalpy::~totalEnthalpy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::totalEnthalpy::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    return true;
}


bool Foam::functionObjects::totalEnthalpy::execute()
{
    const word fieldName(IOobject::groupName("Ha", phaseName_));

    const word thermoName
    (
        IOobject::groupName(fluidThermo::dictName, phaseName_)
    );

    if (mesh_.foundObject<fluidThermo>(thermoName))
    {
        const fluidThermo& thermo = mesh_.lookupObject<fluidThermo>(thermoName);
        const volVectorField& U = mesh_.lookupObject<volVectorField>
        (
            IOobject::groupName("U", phaseName_)
        );

        return store(fieldName, thermo.ha() + 0.5*magSqr(U));
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find fluidThermo " << thermoName
            << " in the database"
            << exit(FatalError);

        return false;
    }
}


bool Foam::functionObjects::totalEnthalpy::write()
{
    return writeLocalObjects::write();
}


// ************************************************************************* //

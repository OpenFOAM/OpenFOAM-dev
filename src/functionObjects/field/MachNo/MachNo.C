/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "MachNo.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(MachNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        MachNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::MachNo::MachNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_),
    phaseName_(word::null),
    UName_("U")
{
    read(dict);
    resetLocalObjectName(IOobject::groupName("Ma", phaseName_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::MachNo::~MachNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::MachNo::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);
    UName_ = dict.lookupOrDefault<word>
    (
        "U",
        IOobject::groupName("U", phaseName_)
    );

    return true;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::MachNo::fields() const
{
    return wordList(UName_);
}


bool Foam::functionObjects::MachNo::execute()
{
    const word fieldName(IOobject::groupName("Ma", phaseName_));

    const word thermoName
    (
        IOobject::groupName(physicalProperties::typeName, phaseName_)
    );

    if (mesh_.foundObject<fluidThermo>(thermoName))
    {
        const fluidThermo& thermo = mesh_.lookupObject<fluidThermo>(thermoName);
        const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

        store(fieldName, mag(U)/sqrt(thermo.gamma()*thermo.p()/thermo.rho()));

        return true;
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


bool Foam::functionObjects::MachNo::write()
{
    return writeLocalObjects::write();
}


// ************************************************************************* //

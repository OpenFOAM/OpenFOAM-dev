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

#include "multiValveEngineState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(multiValveEngineState, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        multiValveEngineState,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMeshMovers::multiValveEngine&
Foam::functionObjects::multiValveEngineState::mve() const
{
    return refCast<const fvMeshMovers::multiValveEngine>(mesh_.mover());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::multiValveEngineState::multiValveEngineState
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::multiValveEngineState::~multiValveEngineState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::multiValveEngineState::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    resetName(typeName);

    return true;
}


void Foam::functionObjects::multiValveEngineState::writeFileHeader
(
    const label i
)
{
    writeHeader(file(), "Engine Motion State");
    writeCommented(file(), "Time");
    writeTabbed(file(), "piston from TDC");
    writeTabbed(file(), "piston speed");
    writeTabbed(file(), "piston displacement");
    writeTabbed(file(), "piston clearance");

    const fvMeshMovers::multiValveEngine& mve = this->mve();
    const fvMeshMovers::multiValveEngine::valveList& valves = mve.valves;

    forAll(valves, valvei)
    {
        const fvMeshMovers::multiValveEngine::valveObject& valve =
            valves[valvei];

        writeTabbed(file(), valve.name + " lift");
        writeTabbed(file(), valve.name + " speed");
        writeTabbed(file(), valve.name + " displacement");
    }

    file() << endl;
}


bool Foam::functionObjects::multiValveEngineState::execute()
{
    return true;
}


bool Foam::functionObjects::multiValveEngineState::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        const fvMeshMovers::multiValveEngine& mve = this->mve();

        const fvMeshMovers::multiValveEngine::pistonObject& piston = mve.piston;

        Log << "Piston: " << nl
            << "    position from TDC: " << piston.position() << " m" << nl
            << "    speed = " << piston.speed() << " m/s" << nl
            << "    displacement: " << piston.displacement() << " m" << nl
            << "    clearance: " << piston.clearance() << " m" << endl;

        writeTime(file());
        file()
            << tab
            << piston.position()  << tab
            << piston.speed() << tab
            << piston.displacement()  << tab
            << piston.clearance();


        const fvMeshMovers::multiValveEngine::valveList& valves = mve.valves;

        forAll(valves, valvei)
        {
            const fvMeshMovers::multiValveEngine::valveObject& valve =
                valves[valvei];

            Log << "Valve " << valve.name << nl
                << "    lift: " << (valve.isOpen() ? valve.lift() : 0)
                << " m "
                << (valve.isOpen() ? "(open)" : "(closed)") << nl
                << "    speed: " << valve.speed() << " m/s" << nl
                << "    displacement: " << valve.displacement() << " m" << endl;

            file()
                << tab
                << (valve.isOpen() ? valve.lift() : 0)  << tab
                << valve.speed() << tab
                << valve.displacement();
        }

        file() << endl;
        Log << endl;
    }

    return true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "systemCall.H"
#include "Time.H"
#include "dynamicCode.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(systemCall, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::systemCall::systemCall
(
    const word& name,
    const objectRegistry&,
    const dictionary& dict,
    const bool
)
:
    name_(name),
    executeCalls_(),
    endCalls_(),
    writeCalls_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::systemCall::~systemCall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::systemCall::read(const dictionary& dict)
{
    dict.readIfPresent("executeCalls", executeCalls_);
    dict.readIfPresent("endCalls", endCalls_);
    dict.readIfPresent("writeCalls", writeCalls_);

    if (executeCalls_.empty() && endCalls_.empty() && writeCalls_.empty())
    {
        WarningIn("Foam::system::read(const dictionary&)")
            << "no executeCalls, endCalls or writeCalls defined."
            << endl;
    }
    else if (!dynamicCode::allowSystemOperations)
    {
        FatalErrorIn("systemCall::read(const dictionary&)")
            << "Executing user-supplied system calls is not enabled by "
            << "default because of " << nl
            << "security issues.  If you trust the case you can enable this "
            << "facility by " << nl
            << "adding to the InfoSwitches setting in the system controlDict:"
            << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system controlDict is either" << nl << nl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/controlDict" << nl << nl
            << "or" << nl << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl << nl
            << exit(FatalError);
    }
}


void Foam::systemCall::execute()
{
    forAll(executeCalls_, callI)
    {
        Foam::system(executeCalls_[callI]);
    }
}


void Foam::systemCall::end()
{
    forAll(endCalls_, callI)
    {
        Foam::system(endCalls_[callI]);
    }
}


void Foam::systemCall::timeSet()
{
    // Do nothing
}


void Foam::systemCall::write()
{
    forAll(writeCalls_, callI)
    {
        Foam::system(writeCalls_[callI]);
    }
}


// ************************************************************************* //

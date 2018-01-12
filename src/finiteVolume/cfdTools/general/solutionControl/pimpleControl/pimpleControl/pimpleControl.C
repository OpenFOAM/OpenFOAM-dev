/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleControl::pimpleControl(fvMesh& mesh, const word& algorithmName)
:
    pimpleNoLoopControl(mesh, algorithmName),
    pimpleLoop(static_cast<solutionControl&>(*this))
{
    read();

    printResidualControls();

    if (nCorrPIMPLE_ > 1)
    {
        printCorrResidualControls(nCorrPIMPLE_);
    }
    else
    {
        Info<< nl << algorithmName << ": Operating solver in PISO mode" << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleControl::~pimpleControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleControl::read()
{
    if (!pimpleNoLoopControl::read() || !pimpleLoop::read())
    {
        return false;
    }

    nCorrPIMPLE_ = dict().lookupOrDefault<label>("nOuterCorrectors", 1);

    return true;
}


bool Foam::pimpleControl::loop()
{
    read();

    if (!pimpleLoop::loop(*this))
    {
        mesh().data::remove("finalIteration");

        return false;
    }

    storePrevIterFields();

    if (finalIter())
    {
        mesh().data::add("finalIteration", true);
    }

    return true;
}


bool Foam::pimpleControl::run(Time& time)
{
    read();

    if (converged())
    {
        time.writeAndEnd();
    }
    else
    {
        storePrevIterFields();
    }

    return time.run();
}


bool Foam::pimpleControl::loop(Time& time)
{
    read();

    if (converged())
    {
        time.writeAndEnd();
    }
    else
    {
        storePrevIterFields();
    }

    return time.loop();
}


// ************************************************************************* //

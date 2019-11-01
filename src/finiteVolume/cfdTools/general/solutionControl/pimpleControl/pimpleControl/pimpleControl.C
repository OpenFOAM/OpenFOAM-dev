/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    pimpleNoLoopControl(mesh, algorithmName, *this),
    pimpleLoop(static_cast<solutionControl&>(*this))
{
    read();

    printResidualControls();

    if (nCorrPimple_ > 1)
    {
        printCorrResidualControls(nCorrPimple_);
    }

    Info<< nl << algorithmName << ": Operating solver in "
        << (mesh.steady() ? "steady-state" : mesh.transient() ? "transient" :
            "mixed steady-state/transient") << " mode with " << nCorrPimple_
        << " outer corrector" << (nCorrPimple_ == 1 ? "" : "s") << nl;

    if (nCorrPimple_ == 1)
    {
        Info<< algorithmName << ": Operating solver in "
            << (mesh.steady() ? "SIMPLE" : "PISO") << " mode" << nl;
    }

    Info<< nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleControl::~pimpleControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleControl::read()
{
    return pimpleNoLoopControl::read() && pimpleLoop::read();
}


bool Foam::pimpleControl::loop()
{
    read();

    if (!pimpleLoop::loop(*this))
    {
        updateFinal();

        return false;
    }

    storePrevIterFields();

    updateFinal();

    return true;
}


bool Foam::pimpleControl::run(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.run();
}


bool Foam::pimpleControl::loop(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.loop();
}


// ************************************************************************* //

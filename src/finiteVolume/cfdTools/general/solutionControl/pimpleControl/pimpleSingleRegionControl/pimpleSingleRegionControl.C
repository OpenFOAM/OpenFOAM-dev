/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "pimpleSingleRegionControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleSingleRegionControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleSingleRegionControl::pimpleSingleRegionControl
(
    pimpleNoLoopControl& pimple
)
:
    pimpleLoop(pimple),
    pimple_(pimple)
{
    pimple_.pimpleLoopPtr_ = this;

    read();

    pimple_.printResidualControls();

    if (nCorr_ > 1)
    {
        pimple_.printCorrResidualControls(nCorr_);
    }

    Info<< nl << pimple_.algorithmName() << ": Operating solver in "
        << (pimple_.mesh().schemes().steady()
        ? "steady-state"
        : pimple_.mesh().schemes().transient() ? "transient" :
            "mixed steady-state/transient") << " mode with " << nCorr_
        << " outer corrector" << (nCorr_ == 1 ? "" : "s") << nl;

    if (nCorr_ == 1)
    {
        Info<< pimple_.algorithmName() << ": Operating solver in "
            << (pimple_.mesh().schemes().steady() ? "SIMPLE" : "PISO")
            << " mode" << nl;
    }

    Info<< nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleSingleRegionControl::~pimpleSingleRegionControl()
{
    pimple_.pimpleLoopPtr_ = nullptr;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleSingleRegionControl::read()
{
    return pimple_.read() && pimpleLoop::read();
}


bool Foam::pimpleSingleRegionControl::loop()
{
    if (!pimpleLoop::loop(pimple_))
    {
        pimple_.updateFinal(pimple_.isFinal(finalIter()));

        return false;
    }

    pimple_.storePrevIterFields();

    pimple_.updateFinal(pimple_.isFinal(finalIter()));

    return true;
}


bool Foam::pimpleSingleRegionControl::run(Time& time)
{
    if (!pimple_.endIfConverged(time))
    {
        pimple_.storePrevIterFields();
    }

    return time.run();
}


bool Foam::pimpleSingleRegionControl::loop(Time& time)
{
    if (!pimple_.endIfConverged(time))
    {
        pimple_.storePrevIterFields();
    }

    return time.loop();
}


// ************************************************************************* //

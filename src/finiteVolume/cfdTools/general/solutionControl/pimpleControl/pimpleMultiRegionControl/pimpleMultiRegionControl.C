/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "pimpleMultiRegionControl.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleMultiRegionControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleMultiRegionControl::pimpleMultiRegionControl
(
    const Time& runTime,
    PtrList<solver>& solvers,
    const word& algorithmName
)
:
    multiRegionSolutionControl(runTime, algorithmName),
    pimpleLoop(static_cast<solutionControl&>(*this)),
    convergenceControl(static_cast<solutionControl&>(*this)),
    correctorConvergenceControl
    (
        static_cast<solutionControl&>(*this),
        "outerCorrector"
    ),
    pimpleControls_(solvers.size()),
    nEcorr_(-1),
    Ecorr_(0)
{
    bool allSteady = true, allTransient = true;

    forAll(solvers, i)
    {
        pimpleControls_.set(i, &solvers[i].pimple);
        pimpleControls_[i].pimpleLoopPtr_ = this;

        allSteady = allSteady && solvers[i].mesh.schemes().steady();
        allTransient = allTransient && solvers[i].mesh.schemes().transient();
    }

    read();

    forAll(solvers, i)
    {
        Info<< nl << algorithmName << ": Region " << solvers[i].mesh.name();
        pimpleControls_[i].printResidualControls();

        if (nCorr_ > 1)
        {
            Info<< nl << algorithmName << ": Region " << solvers[i].mesh.name();
            pimpleControls_[i].printCorrResidualControls(nCorr_);
        }
    }

    Info<< nl << algorithmName << ": Operating solver in "
        << (allSteady ? "steady-state" : allTransient ? "transient" :
            "mixed steady-state/transient") << " mode with " << nCorr_
        << " outer corrector" << (nCorr_ == 1 ? "" : "s") << nl;

    if ((allSteady || allTransient) && nCorr_ == 1)
    {
        Info<< algorithmName << ": Operating solver in "
            << (allSteady ? "SIMPLE" : "PISO") << " mode" << nl;
    }

    Info<< nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleMultiRegionControl::~pimpleMultiRegionControl()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].pimpleLoopPtr_ = nullptr;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleMultiRegionControl::read()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].read();
    }

    nEcorr_ = dict().lookupOrDefault<label>("nEnergyCorrectors", 1);

    if (!pimpleLoop::read())
    {
        return false;
    }

    return true;
}


bool Foam::pimpleMultiRegionControl::hasResidualControls() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        result = result && pimpleControls_[i].hasResidualControls();
    }

    return result;
}


bool Foam::pimpleMultiRegionControl::hasCorrResidualControls() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        result = result && pimpleControls_[i].hasCorrResidualControls();
    }

    return result;
}


bool Foam::pimpleMultiRegionControl::criteriaSatisfied() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        result = pimpleControls_[i].criteriaSatisfied() && result;
    }

    return result;
}


bool Foam::pimpleMultiRegionControl::corrCriteriaSatisfied() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        result = pimpleControls_[i].corrCriteriaSatisfied() && result;
    }

    return result;
}


void Foam::pimpleMultiRegionControl::resetCorrSolveIndex()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].resetCorrSolveIndex();
    }
}


void Foam::pimpleMultiRegionControl::updateCorrSolveIndex()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].updateCorrSolveIndex();
    }
}


bool Foam::pimpleMultiRegionControl::loop()
{
    read();

    if (!pimpleLoop::loop(*this))
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].updateFinal
            (
                pimpleControls_[i].isFinal(finalIter())
            );
        }

        return false;
    }

    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].storePrevIterFields();
    }

    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].updateFinal
        (
            pimpleControls_[i].isFinal(finalIter())
        );
    }

    return true;
}


bool Foam::pimpleMultiRegionControl::correctEnergy()
{
    if (Ecorr_ >= nEcorr_)
    {
        Ecorr_ = 0;
        return false;
    }

    Ecorr_++;

    return true;
}


bool Foam::pimpleMultiRegionControl::run(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].storePrevIterFields();
        }
    }

    return time.run();
}


bool Foam::pimpleMultiRegionControl::loop(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].storePrevIterFields();
        }
    }

    return time.loop();
}


// ************************************************************************* //

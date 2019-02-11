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

#include "pimpleMultiRegionControl.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleMultiRegionControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::Time& Foam::pimpleMultiRegionControl::time
(
    const PtrList<fvMesh>& pimpleMeshes,
    const PtrList<fvMesh>& solidMeshes
)
{
    if (pimpleMeshes.empty() && solidMeshes.empty())
    {
        FatalErrorInFunction
            << "There needs to be at least one region"
            << exit(FatalError);
    }

    if (!pimpleMeshes.empty())
    {
        return pimpleMeshes[0].time();
    }

    return solidMeshes[0].time();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleMultiRegionControl::pimpleMultiRegionControl
(
    PtrList<fvMesh>& pimpleMeshes,
    PtrList<fvMesh>& solidMeshes,
    const word& algorithmName
)
:
    multiRegionSolutionControl(time(pimpleMeshes, solidMeshes), algorithmName),
    pimpleLoop(static_cast<solutionControl&>(*this)),
    convergenceControl(static_cast<solutionControl&>(*this)),
    correctorConvergenceControl
    (
        static_cast<solutionControl&>(*this),
        "outerCorrector"
    ),
    pimpleControls_(),
    solidControls_()
{
    bool allSteady = true, allTransient = true;

    forAll(pimpleMeshes, i)
    {
        pimpleControls_.append
        (
            new pimpleNoLoopControl(pimpleMeshes[i], algorithmName, *this)
        );

        allSteady = allSteady && pimpleMeshes[i].steady();
        allTransient = allTransient && pimpleMeshes[i].transient();
    }

    forAll(solidMeshes, i)
    {
        solidControls_.append
        (
            new solidNoLoopControl(solidMeshes[i], algorithmName, *this)
        );

        allSteady = allSteady && solidMeshes[i].steady();
        allTransient = allTransient && solidMeshes[i].transient();
    }

    read();

    forAll(pimpleMeshes, i)
    {
        Info<< nl << algorithmName << ": Region " << pimpleMeshes[i].name();
        pimpleControls_[i].printResidualControls();

        if (nCorrPimple_ > 1)
        {
            Info<< nl << algorithmName << ": Region " << pimpleMeshes[i].name();
            pimpleControls_[i].printCorrResidualControls(nCorrPimple_);
        }
    }

    forAll(solidMeshes, i)
    {
        Info<< nl << algorithmName << ": Region " << solidMeshes[i].name();
        solidControls_[i].printResidualControls();

        if (nCorrPimple_ > 1)
        {
            Info<< nl << algorithmName << ": Region " << solidMeshes[i].name();
            solidControls_[i].printCorrResidualControls(nCorrPimple_);
        }
    }

    Info<< nl << algorithmName << ": Operating solver in "
        << (allSteady ? "steady-state" : allTransient ? "transient" :
            "mixed steady-state/transient") << " mode with " << nCorrPimple_
        << " outer corrector" << (nCorrPimple_ == 1 ? "" : "s") << nl;

    if ((allSteady || allTransient) && nCorrPimple_ == 1)
    {
        Info<< algorithmName << ": Operating solver in "
            << (allSteady ? "SIMPLE" : "PISO") << " mode" << nl;
    }

    Info<< nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleMultiRegionControl::~pimpleMultiRegionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleMultiRegionControl::read()
{
    forAll(pimpleControls_, i)
    {
        if (!pimpleControls_[i].read())
        {
            return false;
        }
    }
    forAll(solidControls_, i)
    {
        if (!solidControls_[i].read())
        {
            return false;
        }
    }

    const dictionary& solutionDict = dict();

    nCorrPimple_ = solutionDict.lookupOrDefault<label>("nOuterCorrectors", 1);

    return true;
}


bool Foam::pimpleMultiRegionControl::hasResidualControls() const
{
    bool result = true;

    forAll(pimpleControls_, i)
    {
        result = result && pimpleControls_[i].hasResidualControls();
    }
    forAll(solidControls_, i)
    {
        result = result && solidControls_[i].hasResidualControls();
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
    forAll(solidControls_, i)
    {
        result = result && solidControls_[i].hasCorrResidualControls();
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
    forAll(solidControls_, i)
    {
        result = solidControls_[i].criteriaSatisfied() && result;
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
    forAll(solidControls_, i)
    {
        result = solidControls_[i].corrCriteriaSatisfied() && result;
    }

    return result;
}


void Foam::pimpleMultiRegionControl::resetCorrSolveIndex()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].resetCorrSolveIndex();
    }
    forAll(solidControls_, i)
    {
        solidControls_[i].resetCorrSolveIndex();
    }
}


void Foam::pimpleMultiRegionControl::updateCorrSolveIndex()
{
    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].updateCorrSolveIndex();
    }
    forAll(solidControls_, i)
    {
        solidControls_[i].updateCorrSolveIndex();
    }
}


bool Foam::pimpleMultiRegionControl::loop()
{
    read();

    if (!pimpleLoop::loop(*this))
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].updateFinal();
        }
        forAll(solidControls_, i)
        {
            solidControls_[i].updateFinal();
        }

        return false;
    }

    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].storePrevIterFields();
    }
    forAll(solidControls_, i)
    {
        solidControls_[i].storePrevIterFields();
    }

    forAll(pimpleControls_, i)
    {
        pimpleControls_[i].updateFinal();
    }
    forAll(solidControls_, i)
    {
        solidControls_[i].updateFinal();
    }

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
        forAll(solidControls_, i)
        {
            solidControls_[i].storePrevIterFields();
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
        forAll(solidControls_, i)
        {
            solidControls_[i].storePrevIterFields();
        }
    }

    return time.loop();
}


// ************************************************************************* //

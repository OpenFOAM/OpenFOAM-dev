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
    forAll(pimpleMeshes, i)
    {
        pimpleControls_.append
        (
            new pimpleNoLoopControl(pimpleMeshes[i], algorithmName)
        );
    }

    forAll(solidMeshes, i)
    {
        solidControls_.append
        (
            new solidNoLoopControl(solidMeshes[i], algorithmName)
        );
    }

    read();

    forAll(pimpleMeshes, i)
    {
        Info<< nl << algorithmName << ": Region " << pimpleMeshes[i].name();
        pimpleControls_[i].printResidualControls();

        if (nCorrPIMPLE_ > 1)
        {
            Info<< nl << algorithmName << ": Region " << pimpleMeshes[i].name();
            pimpleControls_[i].printCorrResidualControls(nCorrPIMPLE_);
        }
    }

    forAll(solidMeshes, i)
    {
        Info<< nl << algorithmName << ": Region " << solidMeshes[i].name();
        solidControls_[i].printResidualControls();

        if (nCorrPIMPLE_ > 1)
        {
            Info<< nl << algorithmName << ": Region " << solidMeshes[i].name();
            solidControls_[i].printCorrResidualControls(nCorrPIMPLE_);
        }
    }

    if (nCorrPIMPLE_ < 1)
    {
        Info<< nl << algorithmName << ": Operating solver in PISO mode" << nl
            << endl;
    }
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

    nCorrPIMPLE_ = solutionDict.lookupOrDefault<label>("nOuterCorrectors", 1);

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
            pimpleControls_[i].mesh().data::remove("finalIteration");
        }
        forAll(solidControls_, i)
        {
            solidControls_[i].mesh().data::remove("finalIteration");
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

    if (finalIter())
    {
        forAll(pimpleControls_, i)
        {
            pimpleControls_[i].mesh().data::add("finalIteration", true);
        }
        forAll(solidControls_, i)
        {
            solidControls_[i].mesh().data::add("finalIteration", true);
        }
    }

    return true;
}


bool Foam::pimpleMultiRegionControl::run(Time& time)
{
    read();

    if (converged())
    {
        time.writeAndEnd();
    }
    else
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

    if (converged())
    {
        time.writeAndEnd();
    }
    else
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

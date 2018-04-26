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

#include "singleRegionCorrectorConvergenceControl.H"
#include "convergenceControl.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleRegionCorrectorConvergenceControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleRegionCorrectorConvergenceControl::
singleRegionCorrectorConvergenceControl
(
    const singleRegionSolutionControl& control,
    const word& loopName
)
:
    correctorConvergenceControl(control, loopName),
    mesh_(control.mesh()),
    corrResidualControl_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singleRegionCorrectorConvergenceControl::
~singleRegionCorrectorConvergenceControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::singleRegionCorrectorConvergenceControl::readCorrResidualControls()
{
    const dictionary residualDict
    (
        control_.dict().subOrEmptyDict(loopName_ + "ResidualControl")
    );

    DynamicList<corrResidualData> data(corrResidualControl_);

    forAllConstIter(dictionary, residualDict, iter)
    {
        const word& fName = iter().keyword();
        if (!iter().isDict())
        {
            FatalErrorInFunction
                << "Corrector loop convergence criteria specified in "
                << control_.algorithmName() << '.' << residualDict.dictName()
                << " must be given as dictionaries containing \"tolerance\" "
                << "and \"relTol\" entries. Solution convergence criteria are "
                << "specified as single numbers in " << control_.algorithmName()
                << ".residualControl." << exit(FatalError);
        }

        const label fieldi =
            convergenceControl::residualControlIndex
            (
                fName,
                corrResidualControl_,
                false
            );
        if (fieldi == -1)
        {
            corrResidualData rd;
            const dictionary& fieldDict(iter().dict());
            rd.name = fName.c_str();
            rd.absTol = readScalar(fieldDict.lookup("tolerance"));
            rd.relTol = readScalar(fieldDict.lookup("relTol"));
            rd.solveIndex = 0;
            data.append(rd);
        }
        else
        {
            corrResidualData& rd = data[fieldi];
            const dictionary& fieldDict(iter().dict());
            rd.absTol = readScalar(fieldDict.lookup("tolerance"));
            rd.relTol = readScalar(fieldDict.lookup("relTol"));
        }
    }

    corrResidualControl_.transfer(data);

    if (control_.debug > 1)
    {
        forAll(corrResidualControl_, i)
        {
            const corrResidualData& rd = corrResidualControl_[i];
            Info<< residualDict.dictName() << '[' << i << "]:" << nl
                << "    name     : " << rd.name << nl
                << "    absTol   : " << rd.absTol << nl
                << "    relTol   : " << rd.relTol << endl;
        }
    }

    return true;
}


void Foam::singleRegionCorrectorConvergenceControl::printCorrResidualControls
(
    const label n
) const
{
    Info<< nl;

    Info<< control_.algorithmName() << ": "
        << (corrResidualControl_.empty() ? "No c" : "C")
        << "orrector convergence criteria found" << nl;

    forAll(corrResidualControl_, i)
    {
        Info<< control_.algorithmSpace() << "  "
            << corrResidualControl_[i].name << ": tolerance "
            << corrResidualControl_[i].absTol << ", relTol "
            << corrResidualControl_[i].relTol << nl;
    }

    Info<< control_.algorithmSpace() << "  Calclations will do " << n
        << " corrections" << (corrResidualControl_.empty() ? "" :
           " if the convergence criteria are not met") << nl << endl;
}


bool Foam::singleRegionCorrectorConvergenceControl::
hasCorrResidualControls() const
{
    return !corrResidualControl_.empty();
}


bool Foam::singleRegionCorrectorConvergenceControl::
corrCriteriaSatisfied() const
{
    if (!hasCorrResidualControls())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false; // ensure that some checks were actually performed

    if (control_.debug)
    {
        Info<< control_.algorithmName() << ": Correction residuals" << endl;
    }

    const dictionary& solverDict = mesh_.solverPerformanceDict();
    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldi =
            convergenceControl::residualControlIndex
            (
                variableName,
                corrResidualControl_
            );
        if (fieldi != -1)
        {
            scalar firstResidual, residual;
            convergenceControl::getInitialResiduals
            (
                mesh_,
                variableName,
                corrResidualControl_[fieldi].solveIndex,
                iter().stream(),
                firstResidual,
                residual
            );
            const scalar relativeResidual =
                residual/(firstResidual + ROOTVSMALL);

            const bool absCheck =
                residual < corrResidualControl_[fieldi].absTol;
            const bool relCheck =
                relativeResidual < corrResidualControl_[fieldi].relTol;

            checked = true;
            achieved = achieved && (absCheck || relCheck);

            if (control_.debug)
            {
                Info<< control_.algorithmSpace() << "  " << variableName
                    << ": tolerance " << residual << " ("
                    << corrResidualControl_[fieldi].absTol << ")"
                    << ", relTol " << relativeResidual << " ("
                    << corrResidualControl_[fieldi].relTol << ")"
                    << (absCheck || relCheck ? " CONVERGED" : "") << endl;
            }
        }
    }

    return checked && achieved;
}


void Foam::singleRegionCorrectorConvergenceControl::resetCorrSolveIndex()
{
    forAll(corrResidualControl_, i)
    {
        corrResidualControl_[i].solveIndex = 0;
    }
}


void Foam::singleRegionCorrectorConvergenceControl::updateCorrSolveIndex()
{
    const dictionary& solverDict = mesh_.solverPerformanceDict();
    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldi =
            convergenceControl::residualControlIndex
            (
                variableName,
                corrResidualControl_
            );
        if (fieldi != -1)
        {
            getNSolves
            (
                mesh_,
                variableName,
                iter().stream(),
                corrResidualControl_[fieldi].solveIndex
            );
        }
    }
}


// ************************************************************************* //

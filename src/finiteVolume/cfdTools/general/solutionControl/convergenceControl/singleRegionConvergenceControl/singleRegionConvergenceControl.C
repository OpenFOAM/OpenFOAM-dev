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

#include "singleRegionConvergenceControl.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleRegionConvergenceControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleRegionConvergenceControl::singleRegionConvergenceControl
(
    const singleRegionSolutionControl& control
)
:
    convergenceControl(control),
    mesh_(control.mesh()),
    residualControl_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singleRegionConvergenceControl::~singleRegionConvergenceControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::singleRegionConvergenceControl::readResidualControls()
{
    const dictionary residualDict
    (
        control_.dict().subOrEmptyDict("residualControl")
    );

    DynamicList<residualData> data(residualControl_);

    forAllConstIter(dictionary, residualDict, iter)
    {
        const word& fName = iter().keyword();

        if (iter().isDict())
        {
            FatalErrorInFunction
                << "Solution convergence criteria specified in "
                << control_.algorithmName() << '.' << residualDict.dictName()
                << " must be given as single values. Corrector loop "
                << "convergence criteria, if appropriate, are specified as "
                << "dictionaries in " << control_.algorithmName()
                << ".<loopName>ResidualControl." << exit(FatalError);
        }

        const label fieldi =
            residualControlIndex(fName, residualControl_, false);
        if (fieldi == -1)
        {
            residualData rd;
            rd.name = fName.c_str();
            rd.absTol = readScalar(residualDict.lookup(fName));
            data.append(rd);
        }
        else
        {
            residualData& rd = data[fieldi];
            rd.absTol = readScalar(residualDict.lookup(fName));
        }
    }

    residualControl_.transfer(data);

    if (control_.debug > 1)
    {
        forAll(residualControl_, i)
        {
            const residualData& rd = residualControl_[i];
            Info<< residualDict.dictName() << '[' << i << "]:" << nl
                << "    name     : " << rd.name << nl
                << "    absTol   : " << rd.absTol << endl;
        }
    }

    return true;
}


void Foam::singleRegionConvergenceControl::printResidualControls() const
{
    Info<< nl;

    Info<< control_.algorithmName() << ": "
        << (residualControl_.empty() ? "No c" : "C")
        << "onvergence criteria found" << nl;

    forAll(residualControl_, i)
    {
        Info<< control_.algorithmSpace() << "  " << residualControl_[i].name
            << ": tolerance " << residualControl_[i].absTol << nl;
    }

    Info << endl;
}


bool Foam::singleRegionConvergenceControl::hasResidualControls() const
{
    return !residualControl_.empty();
}


bool Foam::singleRegionConvergenceControl::criteriaSatisfied() const
{
    if (!hasResidualControls())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false; // ensure that some checks were actually performed

    if (control_.debug)
    {
        Info<< control_.algorithmName() << ": Residuals" << endl;
    }

    DynamicList<word> fieldNames(getFieldNames(mesh_));

    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];
        const label fieldi =
            residualControlIndex(fieldName, residualControl_);
        if (fieldi != -1)
        {
            scalar residual;
            getInitialResiduals
            (
                mesh_,
                fieldName,
                0,
                residual,
                residual
            );

            checked = true;

            bool absCheck = residual < residualControl_[fieldi].absTol;

            achieved = achieved && absCheck;

            if (control_.debug)
            {
                Info<< control_.algorithmSpace() << "  " << fieldName
                    << ": tolerance " << residual << " ("
                    << residualControl_[fieldi].absTol << ")"
                    << (absCheck ? " CONVERGED" : "") << endl;
            }
        }
    }

    return checked && achieved;
}


// ************************************************************************* //

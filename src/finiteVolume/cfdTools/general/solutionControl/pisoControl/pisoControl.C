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

#include "pisoControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pisoControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pisoControl::pisoControl(fvMesh& mesh, const word& algorithmName)
:
    fluidSolutionControl(mesh, algorithmName),
    nCorrPiso_(-1),
    corrPiso_(0)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pisoControl::~pisoControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pisoControl::read()
{
    if (!fluidSolutionControl::read())
    {
        return false;
    }

    const dictionary& solutionDict = dict();

    nCorrPiso_ = solutionDict.lookupOrDefault<label>("nCorrectors", 1);

    return true;
}


bool Foam::pisoControl::isFinal(const bool finalIter) const
{
    return (finalIter && !anyPisoIter()) || finalPisoIter();
}


bool Foam::pisoControl::correct(const bool finalIter)
{
    read();

    if (finalPisoIter())
    {
        corrPiso_ = 0;

        updateFinal(isFinal(finalIter));

        return false;
    }

    corrPiso_++;

    updateFinal(isFinal(finalIter));

    return true;
}


bool Foam::pisoControl::correctNonOrthogonal(const bool finalIter)
{
    return nonOrthogonalSolutionControl::correctNonOrthogonal
    (
        isFinal(finalIter)
    );
}


bool Foam::pisoControl::run(Time& time)
{
    return time.run();
}


bool Foam::pisoControl::loop(Time& time)
{
    return time.loop();
}


// ************************************************************************* //

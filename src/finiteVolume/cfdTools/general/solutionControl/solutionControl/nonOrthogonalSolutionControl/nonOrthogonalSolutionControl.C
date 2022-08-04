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

#include "nonOrthogonalSolutionControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonOrthogonalSolutionControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonOrthogonalSolutionControl::nonOrthogonalSolutionControl
(
    fvMesh& mesh,
    const word& algorithmName
)
:
    singleRegionSolutionControl(mesh, algorithmName),
    nCorrNonOrth_(-1),
    corrNonOrth_(0)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonOrthogonalSolutionControl::~nonOrthogonalSolutionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::nonOrthogonalSolutionControl::read()
{
    if (!singleRegionSolutionControl::read())
    {
        return false;
    }

    const dictionary& solutionDict = dict();

    nCorrNonOrth_ =
        solutionDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);

    return true;
}


bool Foam::nonOrthogonalSolutionControl::isFinal(const bool finalIter) const
{
    return finalIter && (finalNonOrthogonalIter() || !anyNonOrthogonalIter());
}


bool Foam::nonOrthogonalSolutionControl::correctNonOrthogonal
(
    const bool finalIter
)
{
    read();

    if (finalNonOrthogonalIter())
    {
        corrNonOrth_ = 0;

        updateFinal(isFinal(finalIter));

        return false;
    }

    corrNonOrth_++;

    updateFinal(isFinal(finalIter));

    return true;
}


// ************************************************************************* //

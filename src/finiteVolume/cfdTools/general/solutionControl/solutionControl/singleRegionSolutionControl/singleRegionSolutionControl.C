/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "singleRegionSolutionControl.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleRegionSolutionControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::singleRegionSolutionControl::dependenciesModified() const
{
    return mesh_.solution().modified();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleRegionSolutionControl::singleRegionSolutionControl
(
    fvMesh& mesh,
    const word& algorithmName
)
:
    solutionControl
    (
        mesh,
        mesh.time(),
        (
           !mesh.solution().dict().found(algorithmName)
         && mesh.schemes().steady()
         && mesh.solution().dict().found("SIMPLE")
        )
      ? "SIMPLE"
      : algorithmName
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singleRegionSolutionControl::~singleRegionSolutionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::singleRegionSolutionControl::dict() const
{
    return mesh_.solution().dict().subDict(algorithmName());
}


void Foam::singleRegionSolutionControl::updateFinal
(
    const bool finalIter
) const
{
    mesh_.data::remove("finalIteration");

    if (finalIter)
    {
        mesh_.data::add("finalIteration", true);
    }
}


void Foam::singleRegionSolutionControl::storePrevIterFields()
{
    storePrevIterTypeFields<scalar>();
    storePrevIterTypeFields<vector>();
    storePrevIterTypeFields<sphericalTensor>();
    storePrevIterTypeFields<symmTensor>();
    storePrevIterTypeFields<tensor>();
}


// ************************************************************************* //

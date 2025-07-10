/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2025 OpenFOAM Foundation
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

#include "faceSetSampledSet.H"
#include "faceSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(faceSet, 0);
    addToRunTimeSelectionTable(sampledSet, faceSet, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::faceSet::calcSamples
(
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>&,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    samplingFaces = Foam::faceSet(mesh(), setName_).toc();

    samplingPositions =
        IndirectList<point>(mesh().cellCentres(), samplingFaces);

    samplingSegments = identityMap(samplingFaces.size());

    samplingCells =
        labelList(UIndirectList<label>(mesh().faceOwner(), samplingFaces));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::faceSet::faceSet
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSet(name, mesh, dict),
    setName_(dict.lookup("set"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::faceSet::~faceSet()
{}


// ************************************************************************* //

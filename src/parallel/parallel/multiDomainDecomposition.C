/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "multiDomainDecomposition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiDomainDecomposition, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiDomainDecomposition::multiDomainDecomposition
(
    const processorRunTimes& runTimes,
    const wordList& regionNames
)
:
    multiRegionPrefixer(false, regionNames),
    runTimes_(runTimes),
    regionNames_(regionNames),
    regionMeshes_(regionNames.size())
{
    forAll(regionMeshes_, regioni)
    {
        regionMeshes_.set
        (
            regioni,
            new domainDecomposition(runTimes, regionNames[regioni])
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiDomainDecomposition::~multiDomainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiDomainDecomposition::readDecompose(const bool doSets)
{
    bool result = false;

    forAll(regionMeshes_, regioni)
    {
        if (meshes(regioni)().readDecompose(doSets))
        {
            result = true;

            if (regioni != regionMeshes_.size() - 1)
            {
                Info<< endl;
            }
        }
    }

    return result;
}


bool Foam::multiDomainDecomposition::readReconstruct(const bool doSets)
{
    bool result = false;

    forAll(regionMeshes_, regioni)
    {
        if (meshes(regioni)().readReconstruct(doSets))
        {
            result = true;

            if (regioni != regionMeshes_.size() - 1)
            {
                Info<< endl;
            }
        }
    }

    return result;
}


Foam::fvMesh::readUpdateState
Foam::multiDomainDecomposition::readUpdateDecompose()
{
    fvMesh::readUpdateState result = fvMesh::UNCHANGED;

    forAll(regionMeshes_, regioni)
    {
        const fvMesh::readUpdateState regionResult =
            meshes(regioni)().readUpdateDecompose();

        if
        (
            regioni != regionMeshes_.size() - 1
         && regionResult >= fvMesh::TOPO_CHANGE
        )
        {
            Info<< endl;
        }

        result = result > regionResult ? result : regionResult;
    }

    return result;
}


Foam::fvMesh::readUpdateState
Foam::multiDomainDecomposition::readUpdateReconstruct()
{
    fvMesh::readUpdateState result = fvMesh::UNCHANGED;

    forAll(regionMeshes_, regioni)
    {
        const fvMesh::readUpdateState regionResult =
            meshes(regioni)().readUpdateReconstruct();

        if
        (
            regioni != regionMeshes_.size() - 1
         && regionResult >= fvMesh::TOPO_CHANGE
        )
        {
            Info<< endl;
        }

        result = result > regionResult ? result : regionResult;
    }

    return result;
}


void Foam::multiDomainDecomposition::writeComplete(const bool doSets) const
{
    forAll(regionMeshes_, regioni)
    {
        meshes(regioni)().writeComplete(doSets);
    }
}


void Foam::multiDomainDecomposition::writeProcs(const bool doSets) const
{
    forAll(regionMeshes_, regioni)
    {
        meshes(regioni)().writeProcs(doSets);
    }
}


// ************************************************************************* //

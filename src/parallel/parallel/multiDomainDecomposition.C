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


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::PtrList<Foam::domainDecomposition> Foam::multiDomainDecomposition::init
(
    const processorRunTimes& runTimes,
    const wordList& regionNames
)
{
    Foam::PtrList<domainDecomposition> result(regionNames.size());

    forAll(result, regioni)
    {
        result.set
        (
            regioni,
            new domainDecomposition(runTimes, regionNames[regioni])
        );
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiDomainDecomposition::multiDomainDecomposition
(
    const processorRunTimes& runTimes,
    const wordList& regionNames
)
:
    MultiRegionList<domainDecomposition>(init(runTimes, regionNames)),
    runTimes_(runTimes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiDomainDecomposition::~multiDomainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiDomainDecomposition::readDecompose(const bool doSets)
{
    bool result = false;

    forAll(*this, regioni)
    {
        if (this->operator[](regioni)().readDecompose(doSets))
        {
            result = true;

            if (regioni != size() - 1)
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

    forAll(*this, regioni)
    {
        if (this->operator[](regioni)().readReconstruct(doSets))
        {
            result = true;

            if (regioni != size() - 1)
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

    forAll(*this, regioni)
    {
        const fvMesh::readUpdateState regionResult =
            this->operator[](regioni)().readUpdateDecompose();

        if (regioni != size() - 1 && regionResult >= fvMesh::TOPO_CHANGE)
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

    forAll(*this, regioni)
    {
        const fvMesh::readUpdateState regionResult =
            this->operator[](regioni)().readUpdateReconstruct();

        if (regioni != size() - 1 && regionResult >= fvMesh::TOPO_CHANGE)
        {
            Info<< endl;
        }

        result = result > regionResult ? result : regionResult;
    }

    return result;
}


void Foam::multiDomainDecomposition::writeComplete(const bool doSets) const
{
    forAll(*this, regioni)
    {
        this->operator[](regioni)().writeComplete(doSets);
    }
}


void Foam::multiDomainDecomposition::writeProcs(const bool doSets) const
{
    forAll(*this, regioni)
    {
        this->operator[](regioni)().writeProcs(doSets);
    }
}


// ************************************************************************* //

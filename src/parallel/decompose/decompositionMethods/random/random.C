/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "random.H"
#include "clock.H"
#include "randomGenerator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(random, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        random,
        decomposer
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        random,
        distributor
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::random::random(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict),
    seed_
    (
        decompositionDict.optionalSubDict
        (
            word(decompositionDict.lookup("method")) + "Coeffs"
        ).lookupOrDefault<label>("seed", clock::getTime())
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::random::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    randomGenerator rndGen(seed_);

    // Simple random integer. No guarantees about balance.
    /*
    labelList result(points.size());
    forAll(result, i)
    {
        result[i] = rndGen.sampleAB<label>(0, nDomains());
    }
    */

    // Sorted random scalar. Equal balance.
    scalarList position(points.size());
    forAll(points, i)
    {
        position[i] = rndGen.sample01<scalar>();
    }

    labelList order(points.size());
    sortedOrder(position, order);

    labelList result(points.size());
    forAll(points, i)
    {
        result[i] = (nDomains()*order[i])/points.size();
    }

    return result;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "exponential_oneDimensionalDiscretisation.H"
#include "uniform_oneDimensionalDiscretisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace oneDimensionalDiscretisations
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable
    (
        oneDimensionalDiscretisation,
        exponential,
        dictionary
    );
}
}


// * * * * * * * * * * Private Static Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::oneDimensionalDiscretisations::exponential::coordinates
(
    const dimensionSet& dims,
    const label n,
    const dictionary& dict
)
{
    const scalar min = dict.lookup<scalar>("min", dims);
    const scalar max = dict.lookup<scalar>("max", dims);

    if (min >= max)
    {
        FatalIOErrorInFunction(dict)
            << "The specified minimum and maximum values are not in order"
            << exit(FatalIOError);
    }

    return pow(max/min, uniform::coordinates01(n))*min;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oneDimensionalDiscretisations::exponential::exponential
(
    const word& name,
    const dimensionSet& dims,
    const label n,
    const dictionary& dict
)
:
    oneDimensionalDiscretisation(name, dims, coordinates(dims, n, dict))
{}


// ************************************************************************* //

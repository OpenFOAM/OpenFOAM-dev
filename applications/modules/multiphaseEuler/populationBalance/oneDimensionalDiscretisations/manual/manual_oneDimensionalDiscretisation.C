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

#include "manual_oneDimensionalDiscretisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace oneDimensionalDiscretisations
{
    defineTypeNameAndDebug(manual, 0);
    addToRunTimeSelectionTable
    (
        oneDimensionalDiscretisation,
        manual,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oneDimensionalDiscretisations::manual::manual
(
    const word& name,
    const dimensionSet& dims,
    const label n,
    const dictionary& dict
)
:
    oneDimensionalDiscretisation
    (
        name,
        dims,
        scalarField(dict.lookup<scalarList>("values", dims))
    )
{
    if (coordinates().size() != n)
    {
        FatalIOErrorInFunction(dict)
            << n << " values are required, but " << coordinates().size()
            << " were specified" << exit(FatalIOError);
    }

    for (label i = 0; i < coordinates().size() - 1; ++ i)
    {
        if (coordinates()[i] >= coordinates()[i + 1])
        {
            FatalIOErrorInFunction(dict)
                << "The specified values are not in order"
                << exit(FatalIOError);
        }
    }
}


// ************************************************************************* //

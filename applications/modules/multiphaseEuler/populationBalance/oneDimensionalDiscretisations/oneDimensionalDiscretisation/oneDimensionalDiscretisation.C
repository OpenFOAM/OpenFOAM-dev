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

#include "oneDimensionalDiscretisation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oneDimensionalDiscretisation, 0);
    defineRunTimeSelectionTable(oneDimensionalDiscretisation, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oneDimensionalDiscretisation::oneDimensionalDiscretisation
(
    const word& name,
    const dimensionSet& dims,
    const tmp<scalarField>& coordinates
)
:
    name_(name),
    dims_(dims),
    coordinates_(coordinates)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::oneDimensionalDiscretisation>
Foam::oneDimensionalDiscretisation::New
(
    const word& name,
    const dimensionSet& dims,
    const label n,
    const dictionary& dict
)
{

    const word oddType(dict.lookup("type"));

    Info<< indent << "Selecting one-dimensional discretisation for "
        << name << ": " << oddType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(oddType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown " << typeName << " type "
            << oddType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dims, n, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oneDimensionalDiscretisation::~oneDimensionalDiscretisation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::dimensionedScalar>
Foam::oneDimensionalDiscretisation::dimensionedCoordinates() const
{
    PtrList<dimensionedScalar> result(coordinates().size());

    forAll(coordinates(), i)
    {
        result.set
        (
            i,
            new dimensionedScalar
            (
                name_ + Foam::name(i),
                dims_,
                coordinates()[i]
            )
        );
    }

    return result;
}


// ************************************************************************* //

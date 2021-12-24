/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "fvConstraint.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvConstraint, 0);
    defineRunTimeSelectionTable(fvConstraint, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fvConstraint::constrainType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    return false;
}


template<class Type>
bool Foam::fvConstraint::constrainType(VolField<Type>& field) const
{
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvConstraint::fvConstraint
(
    const word& name,
    const word& constraintType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    constraintType_(constraintType),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(constraintType + "Coeffs"))
{
    Info<< incrIndent << indent << "Name: " << name_
        << endl << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvConstraint> Foam::fvConstraint::New
(
    const word& name,
    const dictionary& coeffs,
    const fvMesh& mesh
)
{
    const word constraintType(coeffs.lookup("type"));

    Info<< indent
        << "Selecting finite volume constraint type " << constraintType << endl;

    libs.open
    (
        coeffs,
        "libs",
        dictionaryConstructorTablePtr_
    );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(constraintType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown fvConstraint " << constraintType << nl << nl
            << "Valid fvConstraints are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<fvConstraint>
    (
        cstrIter()(name, constraintType, coeffs, mesh)
    );
}


Foam::fvConstraint::~fvConstraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fvConstraint::constrainedFields() const
{
    return wordList::null();
}


bool Foam::fvConstraint::constrainsField(const word& fieldName) const
{
    return findIndex(constrainedFields(), fieldName) != -1;
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_CONSTRAINT_CONSTRAIN, fvConstraint);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_CONSTRAINT_CONSTRAIN_FIELD, fvConstraint);


bool Foam::fvConstraint::read(const dictionary& dict)
{
    coeffs_ = dict.optionalSubDict(constraintType_ + "Coeffs");

    return true;
}


// ************************************************************************* //

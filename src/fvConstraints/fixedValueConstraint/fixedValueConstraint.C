/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "fixedValueConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fixedValueConstraint, 0);

    addToRunTimeSelectionTable
    (
        fvConstraint,
        fixedValueConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::fixedValueConstraint::readCoeffs()
{
    fieldValues_.clear();
    forAllConstIter(dictionary, coeffs().subDict("fieldValues"), iter)
    {
        fieldValues_.set
        (
            iter().keyword(),
            new unknownTypeFunction1
            (
                iter().keyword(),
                coeffs().subDict("fieldValues")
            )
        );
    }
}


template<class Type>
bool Foam::fv::fixedValueConstraint::constrainType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const scalar t = mesh().time().value();

    eqn.setValues
    (
        set_.cells(),
        List<Type>(set_.cells().size(), fieldValues_[fieldName]->value<Type>(t))
    );

    return set_.cells().size();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedValueConstraint::fixedValueConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvConstraint(name, modelType, dict, mesh),
    set_(coeffs(), mesh)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::fixedValueConstraint::constrainedFields() const
{
    return fieldValues_.toc();
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_CONSTRAINT_CONSTRAIN,
    fv::fixedValueConstraint
);


void Foam::fv::fixedValueConstraint::updateMesh(const mapPolyMesh& mpm)
{
    set_.updateMesh(mpm);
}


bool Foam::fv::fixedValueConstraint::movePoints()
{
    set_.movePoints();
    return true;
}


bool Foam::fv::fixedValueConstraint::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

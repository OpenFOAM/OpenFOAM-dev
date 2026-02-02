/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "rigidBody.H"
#include "subBody.H"
#include "primitiveFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(rigidBody, 0);
    defineRunTimeSelectionTable(rigidBody, dictionary);
    addToRunTimeSelectionTable(rigidBody, rigidBody, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBD::rigidBody> Foam::RBD::rigidBody::clone() const
{
    return autoPtr<rigidBody>(new rigidBody(*this));
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBD::rigidBody> Foam::RBD::rigidBody::New
(
    const word& name,
    const scalar& m,
    const vector& c,
    const symmTensor& Ic
)
{
    return autoPtr<rigidBody>(new rigidBody(name, m, c, Ic));
}


Foam::autoPtr<Foam::RBD::rigidBody> Foam::RBD::rigidBody::New
(
    const word& name,
    const dictionary& dict
)
{
    const word bodyType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(bodyType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown rigidBody type "
            << bodyType << nl << nl
            << "Valid rigidBody types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<rigidBody>(cstrIter()(name, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBody::~rigidBody()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::RBD::rigidBody::massless() const
{
    return false;
}


void Foam::RBD::rigidBody::merge(const subBody& subBody)
{
    *this = rigidBody
    (
        name(),
        *this + transform(subBody.masterXT(), subBody.body())
    );
}


Foam::tmp<Foam::scalarField> Foam::RBD::rigidBody::sectionMu0s
(
    const direction axis,
    const scalarField& distances
) const
{
    FatalErrorInFunction
        << "Rigid body " << name() << " of type " << type()
        << " does not support the calculation of sectional moments"
        << exit(FatalError);
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::vectorField> Foam::RBD::rigidBody::sectionMu1s
(
    const direction axis,
    const scalarField& distances
) const
{
    sectionMu0s(axis, distances);
    return tmp<vectorField>(nullptr);
}


Foam::tmp<Foam::symmTensorField> Foam::RBD::rigidBody::sectionMu2s
(
    const direction axis,
    const scalarField& distances
) const
{
    sectionMu0s(axis, distances);
    return tmp<symmTensorField>(nullptr);
}


void Foam::RBD::rigidBody::write(Ostream& os) const
{
    writeEntry(os, "type", type());

    writeEntry(os, "mass", m());

    writeEntry(os, "centreOfMass", c());

    writeEntry(os, "inertia", Ic());
}


// ************************************************************************* //

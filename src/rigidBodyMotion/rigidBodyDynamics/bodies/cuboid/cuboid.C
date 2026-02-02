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

#include "cuboid.H"
#include "primitiveFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(cuboid, 0);
    addToRunTimeSelectionTable(rigidBody, cuboid, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBD::rigidBody> Foam::RBD::cuboid::clone() const
{
    return autoPtr<rigidBody>(new cuboid(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::cuboid::~cuboid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::RBD::cuboid::sectionMu0s
(
    const direction axis,
    const scalarField& distances
) const
{
    tmp<scalarField> tResult
    (
        new scalarField(distances.size() - 1, scalar(0))
    );
    scalarField& result = tResult.ref();

    const scalar rho = m()/cmptProduct(L_);
    const scalar cutA = L_[(axis + 1) % 3]*L_[(axis + 2) % 3];

    forAll(result, i)
    {
        const scalar distance0 =
            min(max(distances[i], -L_[axis]/2), +L_[axis]/2);
        const scalar distance1 =
            min(max(distances[i + 1], -L_[axis]/2), +L_[axis]/2);

        result[i] = (distance1 - distance0)*cutA*rho;
    }

    return tResult;
}


Foam::tmp<Foam::vectorField> Foam::RBD::cuboid::sectionMu1s
(
    const direction axis,
    const scalarField& distances
) const
{
    tmp<vectorField> tResult
    (
        new vectorField(distances.size() - 1, vector::zero)
    );
    vectorField& result = tResult.ref();

    const scalar rho = m()/cmptProduct(L_);
    const scalar cutA = L_[(axis + 1) % 3]*L_[(axis + 2) % 3];

    forAll(result, i)
    {
        const scalar distance0 =
            min(max(distances[i], -L_[axis]/2), +L_[axis]/2);
        const scalar distance1 =
            min(max(distances[i + 1], -L_[axis]/2), +L_[axis]/2);

        result[i][axis] = (sqr(distance1) - sqr(distance0))/2*cutA*rho;
    }

    return tResult;
}


Foam::tmp<Foam::symmTensorField> Foam::RBD::cuboid::sectionMu2s
(
    const direction axis,
    const scalarField& distances
) const
{
    tmp<symmTensorField> tResult
    (
        new symmTensorField(distances.size() - 1, symmTensor::zero)
    );
    symmTensorField& result = tResult.ref();

    const scalar rho = m()/cmptProduct(L_);
    const scalar cutA = L_[(axis + 1) % 3]*L_[(axis + 2) % 3];

    static const direction symmTensorDiagComponents[3] =
        { symmTensor::XX, symmTensor::YY, symmTensor::ZZ };

    forAll(result, i)
    {
        const scalar distance0 =
            min(max(distances[i], -L_[axis]/2), +L_[axis]/2);
        const scalar distance1 =
            min(max(distances[i + 1], -L_[axis]/2), +L_[axis]/2);

        result[i][symmTensorDiagComponents[axis]] =
            (pow3(distance1) - pow3(distance0))/3*cutA*rho;

        result[i][symmTensorDiagComponents[(axis + 1) % 3]] =
            pow3(L_[(axis + 1) % 3])/12
           *L_[(axis + 2) % 3]
           *(distance1 - distance0)
           *rho;

        result[i][symmTensorDiagComponents[(axis + 2) % 3]] =
            pow3(L_[(axis + 2) % 3])/12
           *L_[(axis + 1) % 3]
           *(distance1 - distance0)
           *rho;
    }

    return tResult;
}


void Foam::RBD::cuboid::write(Ostream& os) const
{
    writeEntry(os, "type", type());

    writeEntry(os, "mass", m());

    writeEntry(os, "L", L());
}


// ************************************************************************* //

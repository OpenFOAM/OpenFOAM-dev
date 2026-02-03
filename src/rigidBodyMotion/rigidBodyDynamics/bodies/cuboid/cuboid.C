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
    const vector pMin = c() - L_/2;
    const vector pMax = c() + L_/2;

    forAll(result, i)
    {
        vector p0 = pMin, p1 = pMax;
        p0[axis] = min(max(distances[i], pMin[axis]), pMax[axis]);
        p1[axis] = min(max(distances[i + 1], pMin[axis]), pMax[axis]);

        result[i] = cmptProduct(p1 - p0)*rho;
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
    const vector pMin = c() - L_/2;
    const vector pMax = c() + L_/2;

    forAll(result, i)
    {
        vector p0 = pMin, p1 = pMax;
        p0[axis] = min(max(distances[i], pMin[axis]), pMax[axis]);
        p1[axis] = min(max(distances[i + 1], pMin[axis]), pMax[axis]);

        result[i] = cmptProduct(p1 - p0)*rho*(p0 + p1)/2;
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
    const vector pMin = c() - L_/2;
    const vector pMax = c() + L_/2;

    forAll(result, i)
    {
        vector p0 = pMin, p1 = pMax;
        p0[axis] = min(max(distances[i], pMin[axis]), pMax[axis]);
        p1[axis] = min(max(distances[i + 1], pMin[axis]), pMax[axis]);

        const vector d1 = p1 - p0;
        const vector d2 =
        (
            cmptMultiply(p1, p1)
          - cmptMultiply(p0, p0)
        )/2;
        const vector d3 =
        (
            cmptMultiply(p1, cmptMultiply(p1, p1))
          - cmptMultiply(p0, cmptMultiply(p0, p0))
        )/3;

        result[i].xx() = d3.x()*d1.y()*d1.z()*rho;
        result[i].xy() = d2.x()*d2.y()*d1.z()*rho;
        result[i].xz() = d2.x()*d1.y()*d2.z()*rho;
        result[i].yy() = d1.x()*d3.y()*d1.z()*rho;
        result[i].yz() = d1.x()*d2.y()*d2.z()*rho;
        result[i].zz() = d1.x()*d1.y()*d3.z()*rho;
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

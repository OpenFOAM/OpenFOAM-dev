/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
dimensionedTensor dimensionedTensor::T() const
{
    return dimensionedTensor
    (
        name() + ".T()",
        dimensions(),
        value().T()
    );
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

dimensionedScalar tr(const dimensionedTensor& dt)
{
    return dimensionedScalar
    (
        word("tr(", dt.name(), ')'),
        dt.dimensions(),
        tr(dt.value())
    );
}


dimensionedTensor dev(const dimensionedTensor& dt)
{
    return dimensionedTensor
    (
        word("dev(", dt.name(), ')'),
        dt.dimensions(),
        dev(dt.value())
    );
}


dimensionedTensor dev2(const dimensionedTensor& dt)
{
    return dimensionedTensor
    (
        word("dev2(", dt.name(), ')'),
        dt.dimensions(),
        dev2(dt.value())
    );
}


dimensionedScalar det(const dimensionedTensor& dt)
{
    return dimensionedScalar
    (
        word("det(", dt.name(), ')'),
        pow(dt.dimensions(), tensor::dim),
        det(dt.value())
    );
}


dimensionedTensor cof(const dimensionedTensor& dt)
{
    return dimensionedTensor
    (
        word("cof(", dt.name(), ')'),
        pow(dt.dimensions(), tensor::dim - 1),
        cof(dt.value())
    );
}


dimensionedTensor inv(const dimensionedTensor& dt)
{
    return dimensionedTensor
    (
        word("inv(", dt.name(), ')'),
        inv(dt.dimensions()),
        inv(dt.value())
    );
}


dimensionedSymmTensor symm(const dimensionedTensor& dt)
{
    return dimensionedSymmTensor
    (
        word("symm(", dt.name(), ')'),
        dt.dimensions(),
        symm(dt.value())
    );
}

dimensionedSymmTensor twoSymm(const dimensionedTensor& dt)
{
    return dimensionedSymmTensor
    (
        word("twoSymm(", dt.name(), ')'),
        dt.dimensions(),
        twoSymm(dt.value())
    );
}

dimensionedTensor skew(const dimensionedTensor& dt)
{
    return dimensionedTensor
    (
        word("skew(", dt.name(), ')'),
        dt.dimensions(),
        skew(dt.value())
    );
}


dimensionedVector eigenValues(const dimensionedTensor& dt)
{
    return dimensionedVector
    (
        word("eigenValues(", dt.name(), ')'),
        dt.dimensions(),
        eigenValues(dt.value())
    );
}


dimensionedTensor eigenVectors(const dimensionedTensor& dt)
{
    return dimensionedTensor
    (
        word("eigenVectors(", dt.name(), ')'),
        dimless,
        eigenVectors(dt.value())
    );
}


dimensionedVector eigenValues(const dimensionedSymmTensor& dt)
{
    return dimensionedVector
    (
        word("eigenValues(", dt.name(), ')'),
        dt.dimensions(),
        eigenValues(dt.value())
    );
}


dimensionedTensor eigenVectors(const dimensionedSymmTensor& dt)
{
    return dimensionedTensor
    (
        word("eigenVectors(", dt.name(), ')'),
        dimless,
        eigenVectors(dt.value())
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

dimensionedVector operator*(const dimensionedTensor& dt)
{
    return dimensionedVector
    (
        "*"+dt.name(),
        dt.dimensions(),
        *dt.value()
    );
}


dimensionedTensor operator*(const dimensionedVector& dv)
{
    return dimensionedTensor
    (
        "*"+dv.name(),
        dv.dimensions(),
        *dv.value()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

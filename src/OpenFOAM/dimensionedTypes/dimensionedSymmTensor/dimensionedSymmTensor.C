/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "dimensionedSymmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
dimensionedSymmTensor dimensionedSymmTensor::T() const
{
    return dimensionedSymmTensor
    (
        name()+".T()",
        dimensions(),
        value().T()
    );
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

dimensionedSymmTensor sqr(const dimensionedVector& dv)
{
    return dimensionedSymmTensor
    (
        "sqr("+dv.name()+')',
        sqr(dv.dimensions()),
        sqr(dv.value())
    );
}


dimensionedSymmTensor innerSqr(const dimensionedSymmTensor& dt)
{
    return dimensionedSymmTensor
    (
        "innerSqr("+dt.name()+')',
        sqr(dt.dimensions()),
        innerSqr(dt.value())
    );
}


dimensionedScalar tr(const dimensionedSymmTensor& dt)
{
    return dimensionedScalar
    (
        "tr("+dt.name()+')',
        dt.dimensions(),
        tr(dt.value())
    );
}


dimensionedSymmTensor symm(const dimensionedSymmTensor& dt)
{
    return dimensionedSymmTensor
    (
        "symm("+dt.name()+')',
        dt.dimensions(),
        symm(dt.value())
    );
}


dimensionedSymmTensor twoSymm(const dimensionedSymmTensor& dt)
{
    return dimensionedSymmTensor
    (
        "twoSymm("+dt.name()+')',
        dt.dimensions(),
        twoSymm(dt.value())
    );
}


dimensionedSymmTensor dev(const dimensionedSymmTensor& dt)
{
    return dimensionedSymmTensor
    (
        "dev("+dt.name()+')',
        dt.dimensions(),
        dev(dt.value())
    );
}


dimensionedSymmTensor dev2(const dimensionedSymmTensor& dt)
{
    return dimensionedSymmTensor
    (
        "dev2("+dt.name()+')',
        dt.dimensions(),
        dev2(dt.value())
    );
}


dimensionedScalar det(const dimensionedSymmTensor& dt)
{
    return dimensionedScalar
    (
        "det("+dt.name()+')',
        pow(dt.dimensions(), symmTensor::dim),
        det(dt.value())
    );
}


dimensionedSymmTensor cof(const dimensionedSymmTensor& dt)
{
    return dimensionedSymmTensor
    (
        "cof("+dt.name()+')',
        pow(dt.dimensions(), symmTensor::dim - 1),
        cof(dt.value())
    );
}


dimensionedSymmTensor inv(const dimensionedSymmTensor& dt)
{
    return dimensionedSymmTensor
    (
        "inv("+dt.name()+')',
        dimless/dt.dimensions(),
        inv(dt.value())
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

dimensionedVector operator*(const dimensionedSymmTensor& dt)
{
    return dimensionedVector
    (
        "*"+dt.name(),
        dt.dimensions(),
        *dt.value()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

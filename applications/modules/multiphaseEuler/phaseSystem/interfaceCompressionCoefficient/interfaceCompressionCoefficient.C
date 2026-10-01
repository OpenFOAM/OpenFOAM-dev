/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "dimensionedScalar.H"
#include "interfaceCompressionCoefficient.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompressionCoefficient, 0);
    defineBlendedInterfacialModelTypeNameAndDebug
    (
        interfaceCompressionCoefficient,
        0
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompressionCoefficient::~interfaceCompressionCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar&
Foam::interfaceCompressionCoefficient::cAlpha() const
{
    return *this;
}


template<class GeoMesh>
Foam::tmp<Foam::GeometricField<Foam::scalar, GeoMesh>>
Foam::blendedInterfaceCompressionCoefficient::cAlpha() const
{
    return
        evaluate<scalar, GeoMesh>
        (
            &interfaceCompressionCoefficient::cAlpha,
            "cAlpha",
            dimless
        );
}


template Foam::tmp<Foam::volScalarField>
Foam::blendedInterfaceCompressionCoefficient::cAlpha<Foam::fvMesh>() const;


template Foam::tmp<Foam::surfaceScalarField>
Foam::blendedInterfaceCompressionCoefficient::cAlpha<Foam::surfaceMesh>() const;


// ************************************************************************* //

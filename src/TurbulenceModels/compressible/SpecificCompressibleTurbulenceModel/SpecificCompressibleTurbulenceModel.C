/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "SpecificCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicCompressibleTurbulenceModel>
Foam::SpecificCompressibleTurbulenceModel
<
    BasicCompressibleTurbulenceModel
>::SpecificCompressibleTurbulenceModel
(
    const word& type,
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicCompressibleTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicCompressibleTurbulenceModel>
Foam::autoPtr
<
    Foam::SpecificCompressibleTurbulenceModel
    <
        BasicCompressibleTurbulenceModel
    >
>
Foam::SpecificCompressibleTurbulenceModel
<
    BasicCompressibleTurbulenceModel
>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<SpecificCompressibleTurbulenceModel>
    (
        static_cast<SpecificCompressibleTurbulenceModel*>(
        BasicCompressibleTurbulenceModel::New
        (
            geometricOneField(),
            rho,
            U,
            phi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
}


// ************************************************************************* //

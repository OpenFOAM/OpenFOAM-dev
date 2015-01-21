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

#include "SpecificIncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicIncompressibleTurbulenceModel>
Foam::SpecificIncompressibleTurbulenceModel
<
    BasicIncompressibleTurbulenceModel
>::SpecificIncompressibleTurbulenceModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicIncompressibleTurbulenceModel
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

template<class BasicIncompressibleTurbulenceModel>
Foam::autoPtr
<
    Foam::SpecificIncompressibleTurbulenceModel
    <
        BasicIncompressibleTurbulenceModel
    >
>
Foam::SpecificIncompressibleTurbulenceModel
<
    BasicIncompressibleTurbulenceModel
>::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<SpecificIncompressibleTurbulenceModel>
    (
        static_cast<SpecificIncompressibleTurbulenceModel*>(
        BasicIncompressibleTurbulenceModel::New
        (
            geometricOneField(),
            geometricOneField(),
            U,
            phi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
}


// ************************************************************************* //

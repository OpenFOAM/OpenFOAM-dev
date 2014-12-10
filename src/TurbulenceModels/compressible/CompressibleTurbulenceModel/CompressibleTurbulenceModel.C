/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "CompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::CompressibleTurbulenceModel<TransportModel>::
CompressibleTurbulenceModel
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    TurbulenceModel
    <
        geometricOneField,
        volScalarField,
        compressibleTurbulenceModel,
        transportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::CompressibleTurbulenceModel<TransportModel> >
Foam::CompressibleTurbulenceModel<TransportModel>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<CompressibleTurbulenceModel>
    (
        static_cast<CompressibleTurbulenceModel*>(
        TurbulenceModel
        <
            geometricOneField,
            volScalarField,
            compressibleTurbulenceModel,
            transportModel
        >::New
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

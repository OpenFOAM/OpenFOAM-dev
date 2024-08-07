/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMesh.H"

template<class MomentumTransportModel>
inline Foam::autoPtr<MomentumTransportModel> Foam::momentumTransportModel::New
(
    const typename MomentumTransportModel::alphaField& alpha,
    const typename MomentumTransportModel::rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
{
    const word modelType
    (
        IOdictionary
        (
            momentumTransportModel::readModelDict
            (
                U.db(),
                alphaRhoPhi.group()
            )
        ).lookup("simulationType")
    );

    Info<< indent
        << "Selecting turbulence model type " << modelType << endl;

    typename MomentumTransportModel::dictionaryConstructorTable::iterator
        cstrIter =
        MomentumTransportModel::dictionaryConstructorTablePtr_->find(modelType);

    if
    (
        cstrIter
     == MomentumTransportModel::dictionaryConstructorTablePtr_->end()
    )
    {
        FatalErrorInFunction
            << "Unknown " << MomentumTransportModel::typeName << " type "
            << modelType << nl << nl
            << "Valid " << MomentumTransportModel::typeName << " types:" << endl
            << MomentumTransportModel::dictionaryConstructorTablePtr_
               ->sortedToc()
            << exit(FatalError);
    }

    Info<< incrIndent;

    autoPtr<MomentumTransportModel> modelPtr
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, viscosity)
    );

    Info<< decrIndent;

    return modelPtr;
}


// ************************************************************************* //

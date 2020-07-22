/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "parcelCloud.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::parcelCloud> Foam::parcelCloud::New
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g
)
{
    IOdictionary dict
    (
        IOobject
        (
            name + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word type(dict.lookup<word>("type"));

    Info<< "Selecting " << parcelCloud::typeName << " " << type << endl;

    viscosityConstructorTable::iterator cstrIter =
        viscosityConstructorTablePtr_->find(type);

    if (cstrIter == viscosityConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << parcelCloud::typeName << " " << type << nl << nl
            << "Valid cloud types are:" << nl
            << viscosityConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(name, rho, U, mu, g);
}


Foam::autoPtr<Foam::parcelCloud> Foam::parcelCloud::New
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo
)
{
    IOdictionary dict
    (
        IOobject
        (
            name + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word type(dict.lookup<word>("type"));

    Info<< "Selecting " << parcelCloud::typeName << " " << type << endl;

    thermoConstructorTable::iterator cstrIter =
        thermoConstructorTablePtr_->find(type);

    if (cstrIter == thermoConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << parcelCloud::typeName << " " << type << nl << nl
            << "Valid cloud types are:" << nl
            << thermoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(name, rho, U, g, thermo);
}


// ************************************************************************* //

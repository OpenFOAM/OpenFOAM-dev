/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "BirdCarreau.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(BirdCarreau, 0);
    addToRunTimeSelectionTable
    (
        viscosityModel,
        BirdCarreau,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::BirdCarreau::calcNu() const
{
    return
        nuInf_
      + (nu0_ - nuInf_)
       *pow(scalar(1) + pow(k_*strainRate(), a_), (n_ - 1.0)/a_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::BirdCarreau::BirdCarreau
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    BirdCarreauCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    nu0_("nu0", dimViscosity, BirdCarreauCoeffs_),
    nuInf_("nuInf", dimViscosity, BirdCarreauCoeffs_),
    k_("k", dimTime, BirdCarreauCoeffs_),
    n_("n", dimless, BirdCarreauCoeffs_),
    a_(dimensioned<scalar>::lookupOrDefault("a", BirdCarreauCoeffs_, 2)),
    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::BirdCarreau::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    BirdCarreauCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    BirdCarreauCoeffs_.lookup("nu0") >> nu0_;
    BirdCarreauCoeffs_.lookup("nuInf") >> nuInf_;
    BirdCarreauCoeffs_.lookup("k") >> k_;
    BirdCarreauCoeffs_.lookup("n") >> n_;
    a_ = dimensioned<scalar>::lookupOrDefault("a", BirdCarreauCoeffs_, 2);

    return true;
}


// ************************************************************************* //

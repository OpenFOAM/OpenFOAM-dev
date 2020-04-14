/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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

#include "Casson.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Casson, 0);
    addToRunTimeSelectionTable
    (
        viscosityModel,
        Casson,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Casson::calcNu() const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            sqr
            (
                sqrt
                (
                    tau0_
                   /max
                    (
                        strainRate(),
                        dimensionedScalar(dimless/dimTime, vSmall)
                    )
                ) + sqrt(m_)
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::Casson::Casson
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    CassonCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    m_("m", dimViscosity, CassonCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, CassonCoeffs_),
    nuMin_("nuMin", dimViscosity, CassonCoeffs_),
    nuMax_("nuMax", dimViscosity, CassonCoeffs_),
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

bool Foam::viscosityModels::Casson::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    CassonCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    CassonCoeffs_.lookup("m") >> m_;
    CassonCoeffs_.lookup("tau0") >> tau0_;
    CassonCoeffs_.lookup("nuMin_") >> nuMin_;
    CassonCoeffs_.lookup("nuMax_") >> nuMax_;

    return true;
}


// ************************************************************************* //

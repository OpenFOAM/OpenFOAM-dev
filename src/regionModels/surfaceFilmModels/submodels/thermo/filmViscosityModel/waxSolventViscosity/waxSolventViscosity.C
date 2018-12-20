/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "waxSolventViscosity.H"
#include "kinematicSingleLayer.H"
#include "waxSolventEvaporation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(waxSolventViscosity, 0);

addToRunTimeSelectionTable
(
    filmViscosityModel,
    waxSolventViscosity,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void waxSolventViscosity::correctMu()
{
    const kinematicSingleLayer& film = filmType<kinematicSingleLayer>();

    const uniformDimensionedScalarField Wwax
    (
        film.regionMesh().lookupObject<uniformDimensionedScalarField>
        (
            waxSolventEvaporation::typeName + ":Wwax"
        )
    );

    const uniformDimensionedScalarField Wsolvent
    (
        film.regionMesh().lookupObject<uniformDimensionedScalarField>
        (
            waxSolventEvaporation::typeName + ":Wsolvent"
        )
    );

    const uniformDimensionedScalarField Ysolvent0
    (
        film.regionMesh().lookupObject<uniformDimensionedScalarField>
        (
            waxSolventEvaporation::typeName + ":Ysolvent0"
        )
    );

    const volScalarField& Ysolvent
    (
        film.regionMesh().lookupObject<volScalarField>
        (
            waxSolventEvaporation::typeName + ":Ysolvent"
        )
    );

    const volScalarField Xsolvent
    (
        Ysolvent*Wsolvent/((1 - Ysolvent)*Wwax + Ysolvent*Wsolvent)
    );

    const dimensionedScalar Xsolvent0
    (
        Ysolvent0*Wsolvent/((1 - Ysolvent0)*Wwax + Ysolvent0*Wsolvent)
    );

    mu_ = pow(muWax_/muSolvent_, (1 - Xsolvent)/(1 - Xsolvent0))*muSolvent_;
    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waxSolventViscosity::waxSolventViscosity
(
    surfaceFilmRegionModel& film,
    const dictionary& dict,
    volScalarField& mu
)
:
    filmViscosityModel(typeName, film, dict, mu),
    muWax_
    (
        IOobject
        (
            typeName + ":muWax",
            film.regionMesh().time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimDynamicViscosity, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    muWaxModel_
    (
        filmViscosityModel::New
        (
            film,
            coeffDict_.subDict("muWax"),
            muWax_
        )
    ),
    muSolvent_
    (
        IOobject
        (
            typeName + ":muSolvent",
            film.regionMesh().time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimDynamicViscosity, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    muSolventModel_
    (
        filmViscosityModel::New
        (
            film,
            coeffDict_.subDict("muSolvent"),
            muSolvent_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

waxSolventViscosity::~waxSolventViscosity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void waxSolventViscosity::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    muWaxModel_->correct(p, T);
    muSolventModel_->correct(p, T);

    correctMu();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

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

#include "laminar.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "Time.H"
#include "volFields.H"
#include "fvmSup.H"
#include "kinematicSingleLayer.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laminar, 0);
addToRunTimeSelectionTable(momentumTransportModel, laminar, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laminar::laminar
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    momentumTransportModel(type(), film, dict),
    Cf_(coeffDict_.lookup<scalar>("Cf"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

laminar::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volVectorField::Internal> laminar::Us() const
{
    // Evaluate surface velocity assuming parabolic profile
    tmp<volVectorField::Internal> tUs
    (
        volVectorField::Internal::New
        (
            IOobject::modelName("Us", typeName),
            1.5*filmModel_.U()
        )
    );

    return tUs;
}


void laminar::correct()
{}


tmp<fvVectorMatrix> laminar::Su(volVectorField& U) const
{
    // local reference to film model
    const kinematicSingleLayer& film =
        static_cast<const kinematicSingleLayer&>(filmModel_);

    // local references to film fields
    const volScalarField::Internal& mu = film.mu();
    const volScalarField::Internal& rho = film.rho();
    const volScalarField::Internal& delta = film.delta();
    const volVectorField::Internal& Up = film.UPrimary();
    const volScalarField::Internal& rhop = film.rhoPrimary();
    const volScalarField::Internal& VbyA = film.VbyA();

    // Employ simple coeff-based model
    volScalarField::Internal Cs("Cs", Cf_*rhop*mag(Up - U)/VbyA);
    volScalarField::Internal Cw
    (
        "Cw",
        mu/((1.0/3.0)*VbyA*delta + 1e-5*mu*film.time().deltaT()/rho)
    );

    return
    (
       - fvm::Sp(Cs, U) + Cs*Up // surface contribution
       - fvm::Sp(Cw, U) + Cw*film.Uw() // wall contribution
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

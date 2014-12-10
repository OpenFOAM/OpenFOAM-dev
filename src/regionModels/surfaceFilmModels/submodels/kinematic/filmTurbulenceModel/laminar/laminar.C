/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laminar, 0);
addToRunTimeSelectionTable(filmTurbulenceModel, laminar, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laminar::laminar
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    filmTurbulenceModel(type(), owner, dict),
    Cf_(readScalar(coeffDict_.lookup("Cf")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

laminar::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volVectorField> laminar::Us() const
{
    tmp<volVectorField> tUs
    (
        new volVectorField
        (
            IOobject
            (
                typeName + ":Us",
                owner_.regionMesh().time().timeName(),
                owner_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner_.regionMesh(),
            dimensionedVector("zero", dimVelocity, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    // apply quadratic profile
    tUs() = Foam::sqrt(2.0)*owner_.U();
    tUs().correctBoundaryConditions();

    return tUs;
}


tmp<volScalarField> laminar::mut() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":mut",
                owner_.regionMesh().time().timeName(),
                owner_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner_.regionMesh(),
            dimensionedScalar("zero", dimMass/dimLength/dimTime, 0.0)
        )
    );
}


void laminar::correct()
{
    // do nothing
}


tmp<fvVectorMatrix> laminar::Su(volVectorField& U) const
{
    // local reference to film model
    const kinematicSingleLayer& film =
        static_cast<const kinematicSingleLayer&>(owner_);

    // local references to film fields
    const volScalarField& mu = film.mu();
    const volVectorField& Uw = film.Uw();
    const volScalarField& delta = film.delta();
    const volVectorField& Up = film.UPrimary();
    const volScalarField& rhop = film.rhoPrimary();

    // employ simple coeff-based model
    volScalarField Cs("Cs", Cf_*rhop*mag(Up - U));
    volScalarField Cw("Cw", mu/(0.3333*(delta + film.deltaSmall())));
    Cw.min(5000.0);

    return
    (
       - fvm::Sp(Cs, U) + Cs*Up // surface contribution
       - fvm::Sp(Cw, U) + Cw*Uw // wall contribution
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "oneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneEqEddy, 0);
addToRunTimeSelectionTable(LESModel, oneEqEddy, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void oneEqEddy::updateSubGridScaleFields()
{
    nuSgs_ = ck_*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oneEqEddy::oneEqEddy
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    )
{
    bound(k_, kMin_);

    updateSubGridScaleFields();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oneEqEddy::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);

    volScalarField G(GName(), 2.0*nuSgs_*magSqr(symm(gradU)));

    tmp<fvScalarMatrix> kEqn
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(DkEff(), k_)
    ==
       G
     - fvm::Sp(ce_*sqrt(k_)/delta(), k_)
    );

    kEqn().relax();
    kEqn().solve();

    bound(k_, kMin_);

    updateSubGridScaleFields();
}


bool oneEqEddy::read()
{
    if (GenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

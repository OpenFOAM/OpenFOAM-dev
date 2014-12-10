/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "homogeneousDynSmagorinsky.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(homogeneousDynSmagorinsky, 0);
addToRunTimeSelectionTable(LESModel, homogeneousDynSmagorinsky, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void homogeneousDynSmagorinsky::updateSubGridScaleFields
(
    const volSymmTensorField& D
)
{
    nuSgs_ = cD(D)*sqr(delta())*sqrt(magSqr(D));
    nuSgs_.correctBoundaryConditions();
}


dimensionedScalar homogeneousDynSmagorinsky::cD
(
    const volSymmTensorField& D
) const
{
    const volSymmTensorField MM
    (
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D))
    );

    dimensionedScalar MMMM = average(magSqr(MM));

    if (MMMM.value() > VSMALL)
    {
        tmp<volSymmTensorField> LL =
            dev(filter_(sqr(U())) - (sqr(filter_(U()))));

        return 0.5*average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}


dimensionedScalar homogeneousDynSmagorinsky::cI
(
    const volSymmTensorField& D
) const
{
    const volScalarField mm
    (
        sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))))
    );

    dimensionedScalar mmmm = average(magSqr(mm));

    if (mmmm.value() > VSMALL)
    {
        tmp<volScalarField> KK =
            0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

        return average(KK*mm)/mmmm;
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

homogeneousDynSmagorinsky::homogeneousDynSmagorinsky
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

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    bound(k_,  kMin_);

    updateSubGridScaleFields(dev(symm(fvc::grad(U))));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void homogeneousDynSmagorinsky::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    const volSymmTensorField D(dev(symm(gradU)));

    k_ = cI(D)*sqr(delta())*magSqr(D);
    bound(k_,  kMin_);

    updateSubGridScaleFields(D);
}


bool homogeneousDynSmagorinsky::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(coeffDict());

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

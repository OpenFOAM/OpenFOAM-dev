/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "advectionDiffusionPatchDistMethod.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(advectionDiffusion, 0);
    addToRunTimeSelectionTable(patchDistMethod, advectionDiffusion, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::advectionDiffusion::advectionDiffusion
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs),
    coeffs_(dict.optionalSubDict(type() + "Coeffs")),
    pdmPredictor_
    (
        patchDistMethod::New
        (
            coeffs_,
            mesh,
            patchIDs
        )
    ),
    epsilon_(coeffs_.lookupOrDefault<scalar>("epsilon", 0.1)),
    tolerance_(coeffs_.lookupOrDefault<scalar>("tolerance", 1e-3)),
    maxIter_(coeffs_.lookupOrDefault<int>("maxIter", 10)),
    predicted_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::advectionDiffusion::correct(volScalarField& y)
{
    return correct(y, const_cast<volVectorField&>(volVectorField::null()));
}


bool Foam::patchDistMethods::advectionDiffusion::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    if (!predicted_)
    {
        pdmPredictor_->correct(y);
        predicted_ = true;
    }

    volVectorField ny
    (
        IOobject
        (
            "ny",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("ny", dimless, Zero),
        patchTypes<vector>(mesh_, patchIDs_)
    );

    const fvPatchList& patches = mesh_.boundary();
    volVectorField::Boundary& nybf = ny.boundaryFieldRef();

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();
        nybf[patchi] == -patches[patchi].nf();
    }

    int iter = 0;
    scalar initialResidual = 0;

    do
    {
        ny = fvc::grad(y);
        ny /= (mag(ny) + small);

        surfaceVectorField nf(fvc::interpolate(ny));
        nf /= (mag(nf) + small);

        surfaceScalarField yPhi("yPhi", mesh_.Sf() & nf);

        fvScalarMatrix yEqn
        (
            fvm::div(yPhi, y)
          - fvm::Sp(fvc::div(yPhi), y)
          - epsilon_*y*fvm::laplacian(y)
         ==
            dimensionedScalar("1", dimless, 1.0)
        );

        yEqn.relax();
        initialResidual = yEqn.solve().initialResidual();

    } while (initialResidual > tolerance_ && ++iter < maxIter_);

    // Only calculate n if the field is defined
    if (notNull(n))
    {
        n = -ny;
    }

    return true;
}


// ************************************************************************* //

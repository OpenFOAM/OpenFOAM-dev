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

#include "PoissonPatchDistMethod.H"
#include "fvcGrad.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(Poisson, 0);
    addToRunTimeSelectionTable(patchDistMethod, Poisson, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::Poisson::Poisson
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs)
{}


Foam::patchDistMethods::Poisson::Poisson
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::Poisson::correct(volScalarField& y)
{
    return correct(y, const_cast<volVectorField&>(volVectorField::null()));
}


bool Foam::patchDistMethods::Poisson::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    if (!tyPsi_.valid())
    {
        tyPsi_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "yPsi",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("yPsi", sqr(dimLength), 0.0),
                y.boundaryFieldRef().types()
            )
        );
    }
    volScalarField& yPsi = tyPsi_.ref();

    solve(fvm::laplacian(yPsi) == dimensionedScalar("1", dimless, -1.0));

    volVectorField gradyPsi(fvc::grad(yPsi));
    volScalarField magGradyPsi(mag(gradyPsi));

    y = sqrt(magSqr(gradyPsi) + 2*yPsi) - magGradyPsi;

    // Cache yPsi if the mesh is moving otherwise delete
    if (!mesh_.changing())
    {
        tyPsi_.clear();
    }

    // Only calculate n if the field is defined
    if (notNull(n))
    {
        n =
           -gradyPsi
           /max
            (
                magGradyPsi,
                dimensionedScalar("smallMagGradyPsi", dimLength, small)
            );
    }

    return true;
}


// ************************************************************************* //

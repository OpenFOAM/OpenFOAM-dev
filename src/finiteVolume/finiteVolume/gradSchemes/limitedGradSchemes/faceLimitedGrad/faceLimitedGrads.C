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

#include "faceLimitedGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    makeFvGradScheme(faceLimitedGrad)
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::volVectorField>
Foam::fv::faceLimitedGrad<Foam::scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volVectorField> tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    volVectorField& g = tGrad();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    // create limiter
    scalarField limiter(vsf.internalField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar vsfOwn = vsf[own];
        scalar vsfNei = vsf[nei];

        scalar maxFace = max(vsfOwn, vsfNei);
        scalar minFace = min(vsfOwn, vsfNei);
        scalar maxMinFace = rk*(maxFace - minFace);
        maxFace += maxMinFace;
        minFace -= maxMinFace;

        // owner side
        limitFace
        (
            limiter[own],
            maxFace - vsfOwn, minFace - vsfOwn,
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxFace - vsfNei, minFace - vsfNei,
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    const volScalarField::GeometricBoundaryField& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        if (psf.coupled())
        {
            const scalarField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];

                scalar vsfOwn = vsf[own];
                scalar vsfNei = psfNei[pFacei];

                scalar maxFace = max(vsfOwn, vsfNei);
                scalar minFace = min(vsfOwn, vsfNei);
                scalar maxMinFace = rk*(maxFace - minFace);
                maxFace += maxMinFace;
                minFace -= maxMinFace;

                limitFace
                (
                    limiter[own],
                    maxFace - vsfOwn, minFace - vsfOwn,
                    (pCf[pFacei] - C[own]) & g[own]
                );
            }
        }
        else if (psf.fixesValue())
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];

                scalar vsfOwn = vsf[own];
                scalar vsfNei = psf[pFacei];

                scalar maxFace = max(vsfOwn, vsfNei);
                scalar minFace = min(vsfOwn, vsfNei);
                scalar maxMinFace = rk*(maxFace - minFace);
                maxFace += maxMinFace;
                minFace -= maxMinFace;

                limitFace
                (
                    limiter[own],
                    maxFace - vsfOwn, minFace - vsfOwn,
                    (pCf[pFacei] - C[own]) & g[own]
                );
            }
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.internalField() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


template<>
Foam::tmp<Foam::volTensorField>
Foam::fv::faceLimitedGrad<Foam::vector>::calcGrad
(
    const volVectorField& vvf,
    const word& name
) const
{
    const fvMesh& mesh = vvf.mesh();

    tmp<volTensorField> tGrad = basicGradScheme_().calcGrad(vvf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    volTensorField& g = tGrad();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    // create limiter
    scalarField limiter(vvf.internalField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector vvfOwn = vvf[own];
        vector vvfNei = vvf[nei];

        // owner side
        vector gradf = (Cf[facei] - C[own]) & g[own];

        scalar vsfOwn = gradf & vvfOwn;
        scalar vsfNei = gradf & vvfNei;

        scalar maxFace = max(vsfOwn, vsfNei);
        scalar minFace = min(vsfOwn, vsfNei);
        scalar maxMinFace = rk*(maxFace - minFace);
        maxFace += maxMinFace;
        minFace -= maxMinFace;

        limitFace
        (
            limiter[own],
            maxFace - vsfOwn, minFace - vsfOwn,
            magSqr(gradf)
        );


        // neighbour side
        gradf = (Cf[facei] - C[nei]) & g[nei];

        vsfOwn = gradf & vvfOwn;
        vsfNei = gradf & vvfNei;

        maxFace = max(vsfOwn, vsfNei);
        minFace = min(vsfOwn, vsfNei);

        limitFace
        (
            limiter[nei],
            maxFace - vsfNei, minFace - vsfNei,
            magSqr(gradf)
        );
    }


    const volVectorField::GeometricBoundaryField& bvf = vvf.boundaryField();

    forAll(bvf, patchi)
    {
        const fvPatchVectorField& psf = bvf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        if (psf.coupled())
        {
            const vectorField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];

                vector vvfOwn = vvf[own];
                vector vvfNei = psfNei[pFacei];

                vector gradf = (pCf[pFacei] - C[own]) & g[own];

                scalar vsfOwn = gradf & vvfOwn;
                scalar vsfNei = gradf & vvfNei;

                scalar maxFace = max(vsfOwn, vsfNei);
                scalar minFace = min(vsfOwn, vsfNei);
                scalar maxMinFace = rk*(maxFace - minFace);
                maxFace += maxMinFace;
                minFace -= maxMinFace;

                limitFace
                (
                    limiter[own],
                    maxFace - vsfOwn, minFace - vsfOwn,
                    magSqr(gradf)
                );
            }
        }
        else if (psf.fixesValue())
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];

                vector vvfOwn = vvf[own];
                vector vvfNei = psf[pFacei];

                vector gradf = (pCf[pFacei] - C[own]) & g[own];

                scalar vsfOwn = gradf & vvfOwn;
                scalar vsfNei = gradf & vvfNei;

                scalar maxFace = max(vsfOwn, vsfNei);
                scalar minFace = min(vsfOwn, vsfNei);
                scalar maxMinFace = rk*(maxFace - minFace);
                maxFace += maxMinFace;
                minFace -= maxMinFace;

                limitFace
                (
                    limiter[own],
                    maxFace - vsfOwn, minFace - vsfOwn,
                    magSqr(gradf)
                );
            }
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vvf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.internalField() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vvf, g);

    return tGrad;
}


// ************************************************************************* //

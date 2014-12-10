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

#include "cellLimitedGrad.H"
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
    makeFvGradScheme(cellLimitedGrad)
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::volVectorField>
Foam::fv::cellLimitedGrad<Foam::scalar>::calcGrad
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

    scalarField maxVsf(vsf.internalField());
    scalarField minVsf(vsf.internalField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar vsfOwn = vsf[own];
        scalar vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


    const volScalarField::GeometricBoundaryField& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const scalarField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const scalarField maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        //maxVsf *= 1.0/k_;
        //minVsf *= 1.0/k_;
    }


    // create limiter
    scalarField limiter(vsf.internalField().size(), 1.0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                (pCf[pFacei] - C[own]) & g[own]
            );
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
Foam::fv::cellLimitedGrad<Foam::vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volTensorField> tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    volTensorField& g = tGrad();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    vectorField maxVsf(vsf.internalField());
    vectorField minVsf(vsf.internalField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        const vector& vsfOwn = vsf[own];
        const vector& vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


    const volVectorField::GeometricBoundaryField& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchVectorField& psf = bsf[patchi];
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const vectorField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                const vector& vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                const vector& vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const vectorField maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        //maxVsf *= 1.0/k_;
        //minVsf *= 1.0/k_;
    }


    // create limiter
    vectorField limiter(vsf.internalField().size(), vector::one);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                ((pCf[pFacei] - C[own]) & g[own])
            );
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    tensorField& gIf = g.internalField();

    forAll(gIf, celli)
    {
        gIf[celli] = tensor
        (
            cmptMultiply(limiter[celli], gIf[celli].x()),
            cmptMultiply(limiter[celli], gIf[celli].y()),
            cmptMultiply(limiter[celli], gIf[celli].z())
        );
    }

    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


// ************************************************************************* //

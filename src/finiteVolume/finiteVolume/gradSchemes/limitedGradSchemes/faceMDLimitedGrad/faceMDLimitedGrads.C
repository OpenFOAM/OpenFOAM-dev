/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "faceMDLimitedGrad.H"
#include "cellMDLimitedGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvGradScheme(faceMDLimitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::volVectorField>
Foam::fv::faceMDLimitedGrad<Foam::scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volVectorField> tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (k_ < small)
    {
        return tGrad;
    }

    volVectorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar vsfOwn = vsf[own];
        scalar vsfNei = vsf[nei];

        scalar maxFace = max(vsfOwn, vsfNei);
        scalar minFace = min(vsfOwn, vsfNei);

        if (k_ < 1.0)
        {
            scalar maxMinFace = rk*(maxFace - minFace);
            maxFace += maxMinFace;
            minFace -= maxMinFace;
        }

        // owner side
        cellMDLimitedGrad<scalar>::limitFace
        (
            g[own],
            maxFace - vsfOwn,
            minFace - vsfOwn,
            Cf[facei] - C[own]
        );

        // neighbour side
        cellMDLimitedGrad<scalar>::limitFace
        (
            g[nei],
            maxFace - vsfNei,
            minFace - vsfNei,
            Cf[facei] - C[nei]
        );
    }

    const volScalarField::Boundary& bsf = vsf.boundaryField();

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

                if (k_ < 1.0)
                {
                    scalar maxMinFace = rk*(maxFace - minFace);
                    maxFace += maxMinFace;
                    minFace -= maxMinFace;
                }

                cellMDLimitedGrad<scalar>::limitFace
                (
                    g[own],
                    maxFace - vsfOwn,
                    minFace - vsfOwn,
                    pCf[pFacei] - C[own]
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

                if (k_ < 1.0)
                {
                    scalar maxMinFace = rk*(maxFace - minFace);
                    maxFace += maxMinFace;
                    minFace -= maxMinFace;
                }

                cellMDLimitedGrad<scalar>::limitFace
                (
                    g[own],
                    maxFace - vsfOwn,
                    minFace - vsfOwn,
                    pCf[pFacei] - C[own]
                );
            }
        }
    }

    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


template<>
Foam::tmp<Foam::volTensorField>
Foam::fv::faceMDLimitedGrad<Foam::vector>::calcGrad
(
    const volVectorField& vvf,
    const word& name
) const
{
    const fvMesh& mesh = vvf.mesh();

    tmp<volTensorField> tGrad = basicGradScheme_().calcGrad(vvf, name);

    if (k_ < small)
    {
        return tGrad;
    }

    volTensorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector vvfOwn = vvf[own];
        vector vvfNei = vvf[nei];

        vector maxFace = max(vvfOwn, vvfNei);
        vector minFace = min(vvfOwn, vvfNei);

        if (k_ < 1.0)
        {
            vector maxMinFace = rk*(maxFace - minFace);
            maxFace += maxMinFace;
            minFace -= maxMinFace;
        }

        // owner side
        cellMDLimitedGrad<vector>::limitFace
        (
            g[own],
            maxFace - vvfOwn,
            minFace - vvfOwn,
            Cf[facei] - C[own]
        );


        // neighbour side
        cellMDLimitedGrad<vector>::limitFace
        (
            g[nei],
            maxFace - vvfNei,
            minFace - vvfNei,
            Cf[facei] - C[nei]
        );
    }


    const volVectorField::Boundary& bvf = vvf.boundaryField();

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

                vector maxFace = max(vvfOwn, vvfNei);
                vector minFace = min(vvfOwn, vvfNei);

                if (k_ < 1.0)
                {
                    vector maxMinFace = rk*(maxFace - minFace);
                    maxFace += maxMinFace;
                    minFace -= maxMinFace;
                }

                cellMDLimitedGrad<vector>::limitFace
                (
                    g[own],
                    maxFace - vvfOwn, minFace - vvfOwn,
                    pCf[pFacei] - C[own]
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

                vector maxFace = max(vvfOwn, vvfNei);
                vector minFace = min(vvfOwn, vvfNei);

                if (k_ < 1.0)
                {
                    vector maxMinFace = rk*(maxFace - minFace);
                    maxFace += maxMinFace;
                    minFace -= maxMinFace;
                }

                cellMDLimitedGrad<vector>::limitFace
                (
                    g[own],
                    maxFace - vvfOwn,
                    minFace - vvfOwn,
                    pCf[pFacei] - C[own]
                );
            }
        }
    }

    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vvf, g);

    return tGrad;
}


// ************************************************************************* //

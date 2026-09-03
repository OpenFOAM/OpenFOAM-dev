/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "cellMDLimitedGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvGradScheme(cellMDLimitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
void Foam::fv::cellMDLimitedGrad<Foam::scalar>::calcGrad
(
    volInternalVectorField& grad,
    const volScalarField& vsf
) const
{
    basicGradScheme_().calcGrad(grad, vsf);

    if (k_ < small)
    {
        return;
    }

    const fvMesh& mesh = vsf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    scalarField maxVsf(vsf.primitiveField());
    scalarField minVsf(vsf.primitiveField());

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


    const volScalarField::BoundaryField& bsf = vsf.boundaryField();

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

    maxVsf -= vsf.primitiveField();
    minVsf -= vsf.primitiveField();

    if (k_ < 1.0)
    {
        const scalarField maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        // maxVsf *= 1.0/k_;
        // minVsf *= 1.0/k_;
    }


    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            grad[own],
            maxVsf[own],
            minVsf[own],
            Cf[facei] - C[own]
        );

        // neighbour side
        limitFace
        (
            grad[nei],
            maxVsf[nei],
            minVsf[nei],
            Cf[facei] - C[nei]
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
                grad[own],
                maxVsf[own],
                minVsf[own],
                pCf[pFacei] - C[own]
            );
        }
    }
}


template<>
void Foam::fv::cellMDLimitedGrad<Foam::vector>::calcGrad
(
    volInternalTensorField& grad,
    const volVectorField& vvf
) const
{
    basicGradScheme_().calcGrad(grad, vvf);

    if (k_ < small)
    {
        return;
    }

    const fvMesh& mesh = vvf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    vectorField maxVvf(vvf.primitiveField());
    vectorField minVvf(vvf.primitiveField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        const vector& vvfOwn = vvf[own];
        const vector& vvfNei = vvf[nei];

        maxVvf[own] = max(maxVvf[own], vvfNei);
        minVvf[own] = min(minVvf[own], vvfNei);

        maxVvf[nei] = max(maxVvf[nei], vvfOwn);
        minVvf[nei] = min(minVvf[nei], vvfOwn);
    }


    const volVectorField::BoundaryField& bsf = vvf.boundaryField();

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
                const vector& vvfNei = psfNei[pFacei];

                maxVvf[own] = max(maxVvf[own], vvfNei);
                minVvf[own] = min(minVvf[own], vvfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                const vector& vvfNei = psf[pFacei];

                maxVvf[own] = max(maxVvf[own], vvfNei);
                minVvf[own] = min(minVvf[own], vvfNei);
            }
        }
    }

    maxVvf -= vvf.primitiveField();
    minVvf -= vvf.primitiveField();

    if (k_ < 1.0)
    {
        const vectorField maxMinVvf((1.0/k_ - 1.0)*(maxVvf - minVvf));
        maxVvf += maxMinVvf;
        minVvf -= maxMinVvf;

        // maxVvf *= 1.0/k_;
        // minVvf *= 1.0/k_;
    }


    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            grad[own],
            maxVvf[own],
            minVvf[own],
            Cf[facei] - C[own]
        );

        // neighbour side
        limitFace
        (
            grad[nei],
            maxVvf[nei],
            minVvf[nei],
            Cf[facei] - C[nei]
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
                grad[own],
                maxVvf[own],
                minVvf[own],
                pCf[pFacei] - C[own]
            );
        }
    }
}


// ************************************************************************* //

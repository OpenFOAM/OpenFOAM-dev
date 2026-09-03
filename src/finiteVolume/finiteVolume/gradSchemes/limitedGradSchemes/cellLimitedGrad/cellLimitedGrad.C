/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::limitGradient
(
    const Field<scalar>& limiter,
    Field<vector>& gIf
) const
{
    gIf *= limiter;
}


template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::limitGradient
(
    const Field<vector>& limiter,
    Field<tensor>& gIf
) const
{
    forAll(gIf, celli)
    {
        gIf[celli] = tensor
        (
            cmptMultiply(limiter[celli], gIf[celli].x()),
            cmptMultiply(limiter[celli], gIf[celli].y()),
            cmptMultiply(limiter[celli], gIf[celli].z())
        );
    }
}


template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::calcGrad
(
    VolInternalField<typename outerProduct<vector, Type>::type>& grad,
    const VolField<Type>& vf
) const
{
    basicGradScheme_().calcGrad(grad, vf);

    if (k_ < small)
    {
        return;
    }

    const fvMesh& mesh = vf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    Field<Type> maxVf(vf.primitiveField());
    Field<Type> minVf(vf.primitiveField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        const Type& vfOwn = vf[own];
        const Type& vfNei = vf[nei];

        maxVf[own] = max(maxVf[own], vfNei);
        minVf[own] = min(minVf[own], vfNei);

        maxVf[nei] = max(maxVf[nei], vfOwn);
        minVf[nei] = min(minVf[nei], vfOwn);
    }


    const typename VolField<Type>::BoundaryField& bsf =
        vf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchField<Type>& psf = bsf[patchi];
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const Field<Type> psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                const Type& vfNei = psfNei[pFacei];

                maxVf[own] = max(maxVf[own], vfNei);
                minVf[own] = min(minVf[own], vfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                const Type& vfNei = psf[pFacei];

                maxVf[own] = max(maxVf[own], vfNei);
                minVf[own] = min(minVf[own], vfNei);
            }
        }
    }

    maxVf -= vf.primitiveField();
    minVf -= vf.primitiveField();

    if (k_ < 1.0)
    {
        const Field<Type> maxMinVf((1.0/k_ - 1.0)*(maxVf - minVf));
        maxVf += maxMinVf;
        minVf -= maxMinVf;
    }


    // Create limiter initialised to 1
    // Note: the limiter is not permitted to be > 1
    Field<Type> limiter(vf.primitiveField().size(), pTraits<Type>::one);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxVf[own],
            minVf[own],
            (Cf[facei] - C[own]) & grad[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxVf[nei],
            minVf[nei],
            (Cf[facei] - C[nei]) & grad[nei]
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
                maxVf[own],
                minVf[own],
                ((pCf[pFacei] - C[own]) & grad[own])
            );
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    limitGradient(limiter, grad);
}


// ************************************************************************* //

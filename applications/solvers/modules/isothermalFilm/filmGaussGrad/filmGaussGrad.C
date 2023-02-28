/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "filmGaussGrad.H"
#include "filmFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::filmGaussGrad<Type>::calcGrad
(
    const VolField<Type>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<VolField<GradType>> tfGrad
    (
        VolField<GradType>::New
        (
            name,
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    VolField<GradType>& fGrad = tfGrad.ref();

    SurfaceField<Type> ssf(this->tinterpScheme_().interpolate(vsf));

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();
    const scalarField& V = mesh.V();

    Field<GradType>& ifGrad = fGrad;
    const Field<Type>& issf = ssf;

    forAll(owner, facei)
    {
        GradType Sfssf = Sf[facei]*issf[facei];

        ifGrad[owner[facei]] += Sfssf;
        ifGrad[neighbour[facei]] -= Sfssf;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& pFaceCells = p.faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        if (isA<filmFvPatch>(p))
        {
            // Scale the film patch contribution by V/delta

            const scalarField& deltaCoeffs = p.deltaCoeffs();

            forAll(p, facei)
            {
                ifGrad[pFaceCells[facei]] +=
                    0.5*V[pFaceCells[facei]]*deltaCoeffs[facei]
                *(pSf[facei]/mag(pSf[facei]))*pssf[facei];
            }
        }
        else
        {
            forAll(p, facei)
            {
                ifGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
            }
        }
    }

    ifGrad /= mesh.V();

    fGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, fGrad);

    return tfGrad;
}


// ************************************************************************* //

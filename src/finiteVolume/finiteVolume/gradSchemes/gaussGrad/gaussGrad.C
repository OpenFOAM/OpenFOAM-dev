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

#include "gaussGrad.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::gaussGrad<Type>::calcGrad
(
    VolInternalField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type
    >& gGrad,
    const SurfaceField<Type>& sf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = sf.mesh()();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    igGrad = Zero;

    const Field<Type>& isf = sf;

    forAll(owner, facei)
    {
        GradType Sfsf = Sf[facei]*isf[facei];

        igGrad[owner[facei]] += Sfsf;
        igGrad[neighbour[facei]] -= Sfsf;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& pFaceCells = p.faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const fvsPatchField<Type>& psf = sf.boundaryField()[patchi];

        forAll(p, facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*psf[facei];
        }
    }

    igGrad /= mesh.V().primitiveField();
}


template<class Type>
void Foam::fv::gaussGrad<Type>::calcGrad
(
    VolInternalField<typename outerProduct<vector, Type>::type>& grad,
    const VolField<Type>& vf
) const
{
    calcGrad(grad, tinterpScheme_().interpolate(vf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::VolInternalField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type
    >
>
Foam::fv::gaussGrad<Type>::fviGrad
(
    const SurfaceField<Type>& sf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = sf.mesh()();

    tmp<VolInternalField<GradType>> tgGrad
    (
        VolInternalField<GradType>::New
        (
            name,
            mesh,
            sf.dimensions()/dimensions::length
        )
    );

    calcGrad(tgGrad.ref(), sf);

    return tgGrad;
}


template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gaussGrad<Type>::fvcGrad
(
    const SurfaceField<Type>& sf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = sf.mesh()();

    tmp<VolField<GradType>> tgGrad
    (
        VolField<GradType>::New
        (
            name,
            mesh,
            sf.dimensions()/dimensions::length,
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    VolField<GradType>& gGrad = tgGrad.ref();

    calcGrad(tgGrad.ref().internalFieldRef(), sf);
    gGrad.correctBoundaryConditions();

    return tgGrad;
}


// ************************************************************************* //

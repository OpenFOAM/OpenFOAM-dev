/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "fourthGrad.H"
#include "leastSquaresGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::fourthGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    // The fourth-order gradient is calculated in two passes.  First,
    // the standard least-square gradient is assembled.  Then, the
    // fourth-order correction is added to the second-order accurate
    // gradient to complete the accuracy.

    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    // Assemble the second-order least-square gradient
    // Calculate the second-order least-square gradient
    tmp<GeometricField<GradType, fvPatchField, volMesh>> tsecondfGrad
      = leastSquaresGrad<Type>(mesh).grad
        (
            vsf,
            "leastSquaresGrad(" + vsf.name() + ")"
        );
    const GeometricField<GradType, fvPatchField, volMesh>& secondfGrad =
        tsecondfGrad();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tfGrad
    (
        GeometricField<GradType, fvPatchField, volMesh>::New
        (
            name,
            secondfGrad
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& fGrad = tfGrad.ref();

    const vectorField& C = mesh.C();

    const surfaceScalarField& lambda = mesh.weights();

    // Get reference to least square vectors
    const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);
    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    // owner/neighbour addressing
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    // Assemble the fourth-order gradient

    // Internal faces
    forAll(own, facei)
    {
        Type dDotGradDelta = 0.5*
        (
           (C[nei[facei]] - C[own[facei]])
         & (secondfGrad[nei[facei]] - secondfGrad[own[facei]])
        );

        fGrad[own[facei]] -= lambda[facei]*ownLs[facei]*dDotGradDelta;
        fGrad[nei[facei]] -= (1.0 - lambda[facei])*neiLs[facei]*dDotGradDelta;
    }

    // Boundary faces
    forAll(vsf.boundaryField(), patchi)
    {
        if (secondfGrad.boundaryField()[patchi].coupled())
        {
            const fvsPatchVectorField& patchOwnLs =
                ownLs.boundaryField()[patchi];

            const scalarField& lambdap = lambda.boundaryField()[patchi];

            const fvPatch& p = vsf.boundaryField()[patchi].patch();

            const labelUList& faceCells = p.faceCells();

            // Build the d-vectors
            vectorField pd(p.delta());

            const Field<GradType> neighbourSecondfGrad
            (
                secondfGrad.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(faceCells, patchFacei)
            {
                fGrad[faceCells[patchFacei]] -=
                    0.5*lambdap[patchFacei]*patchOwnLs[patchFacei]
                   *(
                        pd[patchFacei]
                      & (
                            neighbourSecondfGrad[patchFacei]
                          - secondfGrad[faceCells[patchFacei]]
                        )
                    );
            }
        }
    }

    fGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, fGrad);

    return tfGrad;
}


// ************************************************************************* //

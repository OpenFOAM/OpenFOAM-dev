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

#include "interpolationCellPointWallModified.H"
#include "syncTools.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class TYPE>
Foam::tmp<Foam::GeometricField<TYPE, Foam::pointPatchField, Foam::pointMesh>>
Foam::interpolationCellPointWallModified<Type>::calcPointField
(
    const GeometricField<TYPE, fvPatchField, volMesh>& psi
) const
{
    FatalErrorInFunction
        << typeName << " interpolation is only defined for vector fields"
        << exit(FatalError);

    return tmp<GeometricField<TYPE, pointPatchField, pointMesh>>(nullptr);
}


template<class Type>
Foam::tmp<Foam::pointVectorField>
Foam::interpolationCellPointWallModified<Type>::calcPointField
(
    const volVectorField& psi
) const
{
    using namespace constant::mathematical;

    const fvMesh& mesh = psi.mesh();


    // Create the point field
    tmp<pointVectorField> tPsip
    (
        pointVectorField::New
        (
            "volPointInterpolateWallModified(" + psi.name() + ')',
            pointMesh::New(mesh),
            dimensioned<Type>("zero", psi.dimensions(), Zero)
        )
    );
    pointVectorField& psip = tPsip.ref();


    // Interpolate to the points with wall patches extrapolated
    {
        wordList patchTypes(psi.boundaryField().size());
        forAll(patchTypes, patchi)
        {
            if (isA<wallPolyPatch>(mesh.boundaryMesh()[patchi]))
            {
                patchTypes[patchi] = zeroGradientFvPatchVectorField::typeName;
            }
            else
            {
                patchTypes[patchi] = calculatedFvPatchVectorField::typeName;
            }
        }
        volVectorField psiExtrapolated
        (
            IOobject
            (
                psi.name() + "Extrapolated",
                mesh.time().timeName(),
                mesh
            ),
            psi,
            patchTypes
        );
        psiExtrapolated.correctBoundaryConditions();
        volPointInterpolation::New(mesh).interpolate(psiExtrapolated, psip);
    }


    // Generate point normals across the entire boundary
    pointField pointNormals(mesh.nPoints(), vector::zero);
    {
        scalarField pointCount(mesh.nPoints(), 0);

        forAll(mesh.boundaryMesh(), patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            forAll(patch, patchFacei)
            {
                const face& f = patch[patchFacei];
                const vector& n = patch.faceNormals()[patchFacei];

                forAll(f, i)
                {
                    pointNormals[f[i]] += n;
                    pointCount[f[i]] += 1;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointNormals,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointCount,
            plusEqOp<scalar>(),
            scalar(0)
        );

        pointNormals /= max(pointCount, small);
    }


    // Calculate the rotation necessary from the vector to the (negative) point
    // normal needed to make the interpolated field point into the mesh
    scalarField theta0(mesh.nPoints(), -vGreat), theta1(mesh.nPoints(), vGreat);
    scalar maxVHatDotN = - vGreat;
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        if (isA<wallPolyPatch>(patch))
        {
            forAll(patch, patchFacei)
            {
                const label facei = patch.start() + patchFacei;
                const face& f = patch[patchFacei];

                for (label i = 1; i < f.size() - 1; ++ i)
                {
                    const triFace tri
                    (
                        tetIndices(mesh.faceOwner()[facei], facei, i)
                       .faceTriIs(mesh)
                    );

                    const vector n = tri.normal(mesh.points());

                    forAll(tri, triPointI)
                    {
                        const label pointi = tri[triPointI];

                        const vector& v = psip[pointi];

                        const scalar vHatDotN = normalised(v) & n;
                        maxVHatDotN = max(maxVHatDotN, vHatDotN);

                        const vector a =
                            normalised(v ^ pointNormals[pointi]);

                        const scalar C = v & n, S = (v ^ a) & n;

                        const scalar theta = atan2(C, - S);

                        theta0[pointi] = max(theta0[pointi], theta);
                        theta1[pointi] = min(theta1[pointi], theta + pi);
                    }
                }
            }
        }
    }

    if (debug)
    {
        reduce(maxVHatDotN, maxOp<scalar>());
        Info<< typeName << ": Maximum in-to-wall dot product before = "
            << maxVHatDotN << endl;
    }

    syncTools::syncPointList
    (
        mesh,
        theta0,
        maxEqOp<scalar>(),
        scalar(0)
    );

    syncTools::syncPointList
    (
        mesh,
        theta1,
        minEqOp<scalar>(),
        scalar(0)
    );


    if (debug > 1)
    {
        pointVectorField(psip.name() + "Before", psip).write();
    }


    // Apply the rotations so that the interpolated field points into the mesh
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        if (isA<wallPolyPatch>(patch))
        {
            forAll(patch.meshPoints(), patchPointi)
            {
                const label pointi = patch.meshPoints()[patchPointi];

                vector& v = psip[pointi];

                if (theta0[pointi] <= 0 && theta1[pointi] >= 0)
                {
                    continue;
                }

                if (theta0[pointi] >= theta1[pointi])
                {
                    v = Zero;
                    continue;
                }

                const scalar theta =
                    theta0[pointi] > 0 ? theta0[pointi] : theta1[pointi];

                const scalar c = cos(theta), s = sin(theta);

                const scalar scale = max(c, 0); // or mag(theta)/(pi/2) ...

                const vector a = normalised(v ^ pointNormals[pointi]);

                v = scale*(tensor::I*c - (*a)*s + sqr(a)*(1 - c)) & v;

                theta0[pointi] = theta1[pointi] = 0;
            }
        }
    }


    // Report the field-normal dot products
    if (debug)
    {
        maxVHatDotN = - vGreat;

        forAll(mesh.boundaryMesh(), patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            if (isA<wallPolyPatch>(patch))
            {
                forAll(patch, patchFacei)
                {
                    const label facei = patch.start() + patchFacei;
                    const face& f = patch[patchFacei];

                    for (label i = 1; i < f.size() - 1; ++ i)
                    {
                        const triFace tri
                        (
                            tetIndices(mesh.faceOwner()[facei], facei, i)
                           .faceTriIs(mesh)
                        );

                        const vector n = tri.normal(mesh.points());

                        forAll(tri, triPointI)
                        {
                            const label pointi = tri[triPointI];

                            const vector& v = psip[pointi];

                            const scalar vHatDotN = normalised(v) & n;
                            maxVHatDotN = max(maxVHatDotN, vHatDotN);
                        }
                    }
                }
            }
        }

        reduce(maxVHatDotN, maxOp<scalar>());
        Info<< typeName << ": Maximum in-to-wall dot product after = "
            << maxVHatDotN << endl;
    }


    // Write the interpolated field
    if (debug > 1)
    {
        psip.write();
    }


    return tPsip;
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationCellPointWallModified<Type>::
interpolationCellPointWallModified
(
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    interpolationCellPoint<Type>(psi, calcPointField(psi))
{}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "volPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "emptyFvPatch.H"
#include "coupledPointPatchField.H"
#include "pointConstraints.H"
#include "UCompactListList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::volPointInterpolation::interpolateUnconstrained
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolateUnconstrained("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const labelListList& pointCells = mesh().pointCells();
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const fvBoundaryMesh& fvbm = mesh().boundary();

    const primitivePatch& boundary = boundaryPtr_();

    // Cache calls to patch coupled flags
    boolList isCoupledPolyPatch(pbm.size(), false);
    boolList isCoupledFvPatch(fvbm.size(), false);
    forAll(isCoupledFvPatch, patchi)
    {
        isCoupledPolyPatch[patchi] = pbm[patchi].coupled();
        isCoupledFvPatch[patchi] = fvbm[patchi].coupled();
    }

    pf = Zero;

    // Interpolate from the cells
    forAll(pointWeights_, pointi)
    {
        forAll(pointWeights_[pointi], pointCelli)
        {
            const label celli = pointCells[pointi][pointCelli];

            pf[pointi] += pointWeights_[pointi][pointCelli]*vf[celli];
        }
    }

    // Get the boundary neighbour field
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary vfBnf
    (
        GeometricField<Type, fvPatchField, volMesh>::null(),
        vf.boundaryField().boundaryNeighbourField()
    );

    // Interpolate from the boundary faces
    forAll(boundaryPointWeights_, bPointi)
    {
        const label pointi = boundary.meshPoints()[bPointi];

        const labelList& pFaces = boundary.pointFaces()[bPointi];

        forAll(boundaryPointWeights_[bPointi], bPointFacei)
        {
            // FV indices
            const labelUList patches =
                mesh().polyBFacePatches()[pFaces[bPointFacei]];
            const labelUList patchFaces =
                mesh().polyBFacePatchFaces()[pFaces[bPointFacei]];

            forAll(boundaryPointWeights_[bPointi][bPointFacei], i)
            {
                // If FV coupled only, add the neighbouring cell value
                if
                (
                    !isCoupledPolyPatch[patches[i]]
                 && isCoupledFvPatch[patches[i]]
                )
                {
                    pf[pointi] +=
                        boundaryPointNbrWeights_[bPointi][bPointFacei][i]
                       *vfBnf[patches[i]][patchFaces[i]];
                }

                // If not coupled, add a weight to the boundary value
                if
                (
                    !isCoupledPolyPatch[patches[i]]
                 && !isCoupledFvPatch[patches[i]]
                )
                {
                    pf[pointi] +=
                        boundaryPointWeights_[bPointi][bPointFacei][i]
                       *vf.boundaryField()[patches[i]][patchFaces[i]];
                }
            }
        }
    }

    // Synchronise
    syncTools::syncPointList(mesh(), pf, plusEqOp<Type>(), pTraits<Type>::zero);
}


template<class Type>
void Foam::volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    interpolateUnconstrained(vf, pf);

    // Apply constraints
    pointConstraints::New(pf.mesh()).constrain(pf);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name,
    const bool cache
) const
{
    typedef GeometricField<Type, pointPatchField, pointMesh> PointFieldType;

    const pointMesh& pm = pointMesh::New(vf.mesh());
    const objectRegistry& db = pm.thisDb();

    if (!cache || vf.mesh().changing())
    {
        // Delete any old occurrences to avoid double registration
        if (db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            PointFieldType& pf =
                db.objectRegistry::template lookupObjectRef<PointFieldType>
                (
                    name
                );

            if (pf.ownedByRegistry())
            {
                solution::cachePrintMessage("Deleting", name, vf);
                pf.release();
                delete &pf;
            }
        }

        tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf
        (
            GeometricField<Type, pointPatchField, pointMesh>::New
            (
                name,
                pm,
                vf.dimensions()
            )
        );

        interpolate(vf, tpf.ref());

        return tpf;
    }
    else
    {
        if (!db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            solution::cachePrintMessage("Calculating and caching", name, vf);
            tmp<PointFieldType> tpf = interpolate(vf, name, false);
            PointFieldType* pfPtr = tpf.ptr();
            regIOobject::store(pfPtr);
            return *pfPtr;
        }
        else
        {
            PointFieldType& pf =
                db.objectRegistry::template lookupObjectRef<PointFieldType>
                (
                    name
                );

            if (pf.upToDate(vf))
            {
                solution::cachePrintMessage("Reusing", name, vf);
                return pf;
            }
            else
            {
                solution::cachePrintMessage("Deleting", name, vf);
                pf.release();
                delete &pf;

                solution::cachePrintMessage("Recalculating", name, vf);
                tmp<PointFieldType> tpf = interpolate(vf, name, false);

                solution::cachePrintMessage("Storing", name, vf);
                PointFieldType* pfPtr = tpf.ptr();
                regIOobject::store(pfPtr);

                return *pfPtr;
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return interpolate(vf, "volPointInterpolate(" + vf.name() + ')', false);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::volPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
) const
{
    tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf =
        interpolate(tvf());

    tvf.clear();

    return tpf;
}


// ************************************************************************* //

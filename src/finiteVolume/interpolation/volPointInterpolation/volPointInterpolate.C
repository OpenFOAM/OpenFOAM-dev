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

#include "volPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "coupledPointPatchField.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::volPointInterpolation::pushUntransformedData
(
    List<Type>& pointData
) const
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh().globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
void Foam::volPointInterpolation::addSeparated
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::addSeparated" << endl;
    }

    typename GeometricField<Type, pointPatchField, pointMesh>::
        Internal& pfi = pf.ref();

    typename GeometricField<Type, pointPatchField, pointMesh>::
        Boundary& pfbf = pf.boundaryFieldRef();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).initSwapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }

    // Block for any outstanding requests
    Pstream::waitRequests();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).swapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }
}


template<class Type>
void Foam::volPointInterpolation::interpolateInternalField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolateInternalField("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const labelListList& pointCells = vf.mesh().pointCells();

    // Multiply volField by weighting factor matrix to create pointField
    forAll(pointCells, pointi)
    {
        if (!isPatchPoint_[pointi])
        {
            const scalarList& pw = pointWeights_[pointi];
            const labelList& ppc = pointCells[pointi];

            pf[pointi] = Zero;

            forAll(ppc, pointCelli)
            {
                pf[pointi] += pw[pointCelli]*vf[ppc[pointCelli]];
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volPointInterpolation::flatBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const fvBoundaryMesh& bm = mesh.boundary();

    tmp<Field<Type>> tboundaryVals
    (
        new Field<Type>(mesh.nFaces()-mesh.nInternalFaces())
    );
    Field<Type>& boundaryVals = tboundaryVals.ref();

    forAll(vf.boundaryField(), patchi)
    {
        label bFacei = bm[patchi].patch().start() - mesh.nInternalFaces();

        if
        (
           !isA<emptyFvPatch>(bm[patchi])
        && !vf.boundaryField()[patchi].coupled()
        )
        {
            SubList<Type>
            (
                boundaryVals,
                vf.boundaryField()[patchi].size(),
                bFacei
            ) = vf.boundaryField()[patchi];
        }
        else
        {
            const polyPatch& pp = bm[patchi].patch();

            forAll(pp, i)
            {
                boundaryVals[bFacei++] = Zero;
            }
        }
    }

    return tboundaryVals;
}


template<class Type>
void Foam::volPointInterpolation::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    const primitivePatch& boundary = boundaryPtr_();

    Field<Type>& pfi = pf.primitiveFieldRef();

    // Get face data in flat list
    tmp<Field<Type>> tboundaryVals(flatBoundaryField(vf));
    const Field<Type>& boundaryVals = tboundaryVals();

    // Do points on 'normal' patches from the surrounding patch faces
    forAll(boundary.meshPoints(), i)
    {
        const label pointi = boundary.meshPoints()[i];

        if (isPatchPoint_[pointi])
        {
            const labelList& pFaces = boundary.pointFaces()[i];
            const scalarList& pWeights = boundaryPointWeights_[i];

            Type& val = pfi[pointi];

            val = Zero;
            forAll(pFaces, j)
            {
                if (boundaryIsPatchFace_[pFaces[j]])
                {
                    val += pWeights[j]*boundaryVals[pFaces[j]];
                }
            }
        }
    }

    // Sum collocated contributions
    pointConstraints::syncUntransformedData(mesh(), pfi, plusEqOp<Type>());

    // And add separated contributions
    addSeparated(pf);

    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves.
    pushUntransformedData(pfi);


    // Detect whether the field has overridden constraint patch types. If not,
    // we are done, so return.
    bool havePatchTypes = false;
    wordList patchTypes(pf.boundaryField().size(), word::null);
    forAll(pf.boundaryField(), patchi)
    {
        const word patchType = pf.boundaryField()[patchi].patchType();
        if (patchType != word::null)
        {
            havePatchTypes = true;
            patchTypes[patchi] = patchType;
        }
    }
    if (!havePatchTypes)
    {
        return;
    }

    // If the patch types have been overridden than we need to re-normalise the
    // boundary points weights. Re-calculate the weight sum.
    pointScalarField psw
    (
        IOobject
        (
            "volPointSumWeights",
            mesh().polyMesh::instance(),
            mesh()
        ),
        pointMesh::New(mesh()),
        dimensionedScalar(dimless, 0),
        wordList
        (
            pf.boundaryField().size(),
            calculatedPointPatchField<scalar>::typeName
        ),
        patchTypes
    );

    scalarField& pswi = psw.primitiveFieldRef();

    forAll(boundary.meshPoints(), i)
    {
        const label pointi = boundary.meshPoints()[i];

        if (isPatchPoint_[pointi])
        {
            pswi[pointi] = sum(boundaryPointWeights_[i]);
        }
    }

    pointConstraints::syncUntransformedData(mesh(), pswi, plusEqOp<scalar>());

    addSeparated(psw);

    pushUntransformedData(pswi);

    // Apply the new weight sum to the result to re-normalise it
    forAll(boundary.meshPoints(), i)
    {
        const label pointi = boundary.meshPoints()[i];

        if (isPatchPoint_[pointi])
        {
            pfi[pointi] /= pswi[pointi];
        }
    }
}


template<class Type>
void Foam::volPointInterpolation::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const bool overrideFixedValue
) const
{
    interpolateBoundaryField(vf, pf);

    // Apply constraints
    const pointConstraints& pcs = pointConstraints::New(pf.mesh());

    pcs.constrain(pf, overrideFixedValue);
}


template<class Type>
void Foam::volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    interpolateInternalField(vf, pf);

    // Interpolate to the patches preserving fixed value BCs
    interpolateBoundaryField(vf, pf, false);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const wordList& patchFieldTypes
) const
{
    const pointMesh& pm = pointMesh::New(vf.mesh());

    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf
    (
        GeometricField<Type, pointPatchField, pointMesh>::New
        (
            "volPointInterpolate(" + vf.name() + ')',
            pm,
            vf.dimensions(),
            patchFieldTypes
        )
    );

    interpolateInternalField(vf, tpf());

    // Interpolate to the patches overriding fixed value BCs
    interpolateBoundaryField(vf, tpf(), true);

    return tpf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::volPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const wordList& patchFieldTypes
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf =
        interpolate(tvf(), patchFieldTypes);
    tvf.clear();
    return tpf;
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

            if (pf.upToDate(vf))    // TBD: , vf.mesh().points()))
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

                // Note: return reference, not pointer
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
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf =
        interpolate(tvf());
    tvf.clear();
    return tpf;
}


// ************************************************************************* //

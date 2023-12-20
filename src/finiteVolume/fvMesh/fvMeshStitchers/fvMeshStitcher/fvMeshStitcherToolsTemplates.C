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

#include "fvMeshStitcherTools.H"
#include "surfaceFields.H"
#include "coupledFvPatch.H"
#include "nonConformalBoundary.H"
#include "nonConformalFvPatch.H"
#include "nonConformalErrorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvMeshStitcherTools::fieldMap
(
    const Field<Type>& f,
    const labelUList& addr,
    const label addr0
)
{
    tmp<Field<Type>> tresult(new Field<Type>(addr.size()));
    forAll(addr, i)
    {
        tresult.ref()[i] = f[addr[i] - addr0];
    }
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvMeshStitcherTools::fieldMap
(
    const tmp<Field<Type>>& f,
    const labelUList& addr,
    const label addr0
)
{

    tmp<Field<Type>> tresult = fieldMap(f(), addr, addr0);
    f.clear();
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvMeshStitcherTools::fieldRMapSum
(
    const Field<Type>& f,
    const label size,
    const labelUList& addr,
    const label addr0
)
{
    tmp<Field<Type>> tresult(new Field<Type>(size, Zero));
    forAll(addr, i)
    {
        tresult.ref()[addr[i] - addr0] += f[i];
    }
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvMeshStitcherTools::fieldRMapSum
(
    const tmp<Field<Type>>& f,
    const label size,
    const labelUList& addr,
    const label addr0
)
{
    tmp<Field<Type>> tresult = fieldRMapSum(f(), size, addr, addr0);
    f.clear();
    return tresult;
}


template<class Type>
Foam::tmp<Foam::fvMeshStitcherTools::SurfaceFieldBoundary<Type>>
Foam::fvMeshStitcherTools::conformedNcBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb
)
{
    const bool isFluxField = isFlux(fieldb[0].internalField());

    const fvBoundaryMesh& fvbm = fieldb[0].patch().boundaryMesh();

    tmp<SurfaceFieldBoundary<Type>> tncFieldb
    (
        new SurfaceFieldBoundary<Type>
        (
            SurfaceField<Type>::Internal::null(),
            fieldb
        )
    );
    SurfaceFieldBoundary<Type>& ncFieldb = tncFieldb.ref();

    ncFieldb == pTraits<Type>::zero;

    const surfaceScalarField::Boundary origNcMagSfb
    (
        surfaceScalarField::Internal::null(),
        fvMeshStitcherTools::origNcMagSfb(fvbm.mesh())
    );

    const labelList origPatchIDs =
        nonConformalBoundary::New(fvbm.mesh()).allOrigPatchIndices();

    // Accumulate the non-conformal parts of the field into the original faces
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvbm[ncPatchi]);

        const label origPatchi = ncFvp.origPatchIndex();
        const fvPatch& origFvp = ncFvp.origPatch();

        const scalarField& ncNcMagSf = ncFvp.patch().magSf();

        // Sum properties with an area-weight, unless this is a flux. Fluxes
        // already scale with the area.
        ncFieldb[origPatchi] +=
            fvMeshStitcherTools::fieldRMapSum
            (
                isFluxField ? fieldb[ncPatchi] : ncNcMagSf*fieldb[ncPatchi],
                origFvp.size(),
                ncFvp.polyFaces(),
                origFvp.start()
            );
    }

    // Scale or average as appropriate
    forAll(origPatchIDs, i)
    {
        const label origPatchi = origPatchIDs[i];
        const fvPatch& origFvp = fvbm[origPatchi];

        // If this is a flux then scale to the total size of the face
        if (isFluxField)
        {
            const scalarField origTotalMagSf
            (
                origFvp.magSf() + origNcMagSfb[origPatchi]
            );

            ncFieldb[origPatchi] *=
                origTotalMagSf
               /max(origNcMagSfb[origPatchi], small*origTotalMagSf);
        }
        // If this is not a flux then divide by the area to create an
        // area-weighted average
        else
        {
            ncFieldb[origPatchi] /=
                max(origNcMagSfb[origPatchi], vSmall);
        }
    }

    return tncFieldb;
}


template<class Type>
Foam::tmp<Foam::fvMeshStitcherTools::SurfaceFieldBoundary<Type>>
Foam::fvMeshStitcherTools::conformedOrigBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb
)
{
    const bool isFluxField = isFlux(fieldb[0].internalField());

    const fvBoundaryMesh& fvbm = fieldb[0].patch().boundaryMesh();

    tmp<SurfaceFieldBoundary<Type>> torigFieldb
    (
        new SurfaceFieldBoundary<Type>
        (
            SurfaceField<Type>::Internal::null(),
            fieldb
        )
    );
    SurfaceFieldBoundary<Type>& origFieldb = torigFieldb.ref();

    const surfaceScalarField::Boundary origNcMagSfb
    (
        surfaceScalarField::Internal::null(),
        fvMeshStitcherTools::origNcMagSfb(fvbm.mesh())
    );

    const labelList origPatchIDs =
        nonConformalBoundary::New(fvbm.mesh()).allOrigPatchIndices();

    // Scale or average as appropriate
    forAll(origPatchIDs, i)
    {
        const label origPatchi = origPatchIDs[i];
        const fvPatch& origFvp = fvbm[origPatchi];

        // If this is a flux then scale to the total size of the face
        if (isFluxField)
        {
            const scalarField origTotalMagSf
            (
                origFvp.magSf() + origNcMagSfb[origPatchi]
            );

            origFieldb[origPatchi] *=
                origTotalMagSf
               /max(origFvp.magSf(), small*origTotalMagSf);
        }
    }

    return torigFieldb;
}


template<class Type>
Foam::tmp<Foam::fvMeshStitcherTools::SurfaceFieldBoundary<Type>>
Foam::fvMeshStitcherTools::unconformedBoundaryField
(
    const SurfaceFieldBoundary<Type>& ncFieldb,
    const SurfaceFieldBoundary<Type>& origFieldb
)
{
    const bool isFluxField = isFlux(origFieldb[0].internalField());

    const fvBoundaryMesh& fvbm = origFieldb[0].patch().boundaryMesh();

    const surfaceScalarField::Boundary& magSfb =
        fvbm.mesh().magSf().boundaryField();

    // Initialise the result and copy the original fields
    tmp<SurfaceFieldBoundary<Type>> tfieldb
    (
        new SurfaceFieldBoundary<Type>
        (
            SurfaceField<Type>::Internal::null(),
            origFieldb
        )
    );
    SurfaceFieldBoundary<Type>& fieldb = tfieldb.ref();

    // Map the conformed non-conformal values into the non-conformal patch
    // fields
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        const label origPatchi = ncFvp.origPatchIndex();
        const fvPatch& origFvp = ncFvp.origPatch();

        fieldb[ncPatchi] =
            fvMeshStitcherTools::fieldMap
            (
                ncFieldb[origPatchi],
                ncFvp.polyFaces(),
                origFvp.start()
            );
    }

    const labelList origPatchIDs =
        nonConformalBoundary::New(fvbm.mesh()).allOrigPatchIndices();

    // If a flux then scale down to the part face area
    if (isFluxField)
    {
        const surfaceScalarField::Boundary origNcMagSfb
        (
            surfaceScalarField::Internal::null(),
            fvMeshStitcherTools::origNcMagSfb(fvbm.mesh())
        );

        // Scale the original patch fields
        forAll(origPatchIDs, i)
        {
            const label origPatchi = origPatchIDs[i];

            const scalarField origTotalMagSf
            (
                magSfb[origPatchi] + origNcMagSfb[origPatchi]
            );

            fieldb[origPatchi] *= magSfb[origPatchi]/origTotalMagSf;
        }

        // Scale the non-conformal patch fields
        forAll(fvbm, ncPatchi)
        {
            const fvPatch& fvp = fvbm[ncPatchi];

            if (!isA<nonConformalFvPatch>(fvp)) continue;

            const nonConformalFvPatch& ncFvp =
                refCast<const nonConformalFvPatch>(fvp);

            const label origPatchi = ncFvp.origPatchIndex();
            const fvPatch& origFvp = ncFvp.origPatch();

            const scalarField ncTotalMagSf
            (
                fvMeshStitcherTools::fieldMap
                (
                    magSfb[origPatchi] + origNcMagSfb[origPatchi],
                    ncFvp.polyFaces(),
                    origFvp.start()
                )
            );

            fieldb[ncPatchi] *= magSfb[ncPatchi]/ncTotalMagSf;
        }
    }

    // Overwrite error values
    forAll(fvbm, errorPatchi)
    {
        const fvPatch& fvp = fvbm[errorPatchi];

        if (!isA<nonConformalErrorFvPatch>(fvp)) continue;

        const nonConformalErrorFvPatch& errorFvp =
            refCast<const nonConformalErrorFvPatch>(fvp);

        const label origPatchi = errorFvp.origPatchIndex();
        const fvPatch& origFvp = errorFvp.origPatch();

        if (isFluxField)
        {
            fieldb[errorPatchi] = Zero;
        }
        else
        {
            fieldb[errorPatchi] =
                fvMeshStitcherTools::fieldMap
                (
                    origFieldb[origPatchi],
                    errorFvp.polyFaces(),
                    origFvp.start()
                );
        }
    }

    return tfieldb;
}


template<class Type>
Foam::tmp<Foam::fvMeshStitcherTools::SurfaceFieldBoundary<Type>>
Foam::fvMeshStitcherTools::synchronisedBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb,
    const bool flip,
    const scalar ownerWeight,
    const scalar neighbourWeight
)
{
    const fvBoundaryMesh& fvbm = fieldb[0].patch().boundaryMesh();

    tmp<SurfaceFieldBoundary<Type>> tsyncFieldb
    (
        new SurfaceFieldBoundary<Type>
        (
            SurfaceField<Type>::Internal::null(),
            fieldb
        )
    );
    SurfaceFieldBoundary<Type>& syncFieldb = tsyncFieldb.ref();

    SurfaceFieldBoundary<Type> fieldbNbr
    (
        SurfaceField<Type>::Internal::null(),
        fieldb.boundaryNeighbourField()
    );

    forAll(fvbm, patchi)
    {
        const fvPatch& fvp = fvbm[patchi];

        if (!fvp.coupled()) continue;

        const coupledFvPatch& cFvp = refCast<const coupledFvPatch>(fvp);

        const scalar w = cFvp.owner() ? ownerWeight : neighbourWeight;
        const scalar v = cFvp.owner() ? neighbourWeight : ownerWeight;

        syncFieldb[patchi] =
            w*syncFieldb[patchi] + (flip ? -v : +v)*fieldbNbr[patchi];
    }

    return tsyncFieldb;
}


template<class Type>
Foam::tmp<Foam::fvMeshStitcherTools::SurfaceFieldBoundary<Type>>
Foam::fvMeshStitcherTools::synchronisedBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb
)
{
    const bool isFluxField = isFlux(fieldb[0].internalField());

    return synchronisedBoundaryField<Type>
    (
        fieldb,
        isFluxField,
        0.5,
        0.5
    );
}


// ************************************************************************* //

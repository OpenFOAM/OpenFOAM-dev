/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "conformedFvsPatchField.H"
#include "fvMeshStitcherTools.H"
#include "nonConformalBoundary.H"
#include "nonConformalFvPatch.H"
#include "nonConformalErrorFvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * * //

template<class Type>
Foam::conformedFvsPatchField<Type>::conformedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    autoPtr<fvsPatchField<Type>>&& origFieldPtr,
    autoPtr<calculatedFvsPatchField<Type>>&& ncFieldPtr
)
:
    fvsPatchField<Type>(p, iF),
    origFieldPtr_(origFieldPtr),
    ncFieldPtr_(ncFieldPtr)
{}


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
void Foam::conformedFvsPatchField<Type>::conform
(
    typename SurfaceField<Type>::Boundary& bF
)
{
    const DimensionedField<Type, surfaceMesh>& iF = bF[0].internalField();

    const fvBoundaryMesh& fvbm = iF.mesh().boundary();

    const labelList origPatchIndices =
        nonConformalBoundary::New(iF.mesh()).allOrigPatchIndices();

    // Evaluate the conformed orig and non-conformal boundary fields
    const typename SurfaceField<Type>::Boundary origBf
    (
        SurfaceField<Type>::Internal::null(),
        fvMeshStitcherTools::conformedOrigBoundaryField(bF)
    );
    const typename SurfaceField<Type>::Boundary ncBf
    (
        SurfaceField<Type>::Internal::null(),
        fvMeshStitcherTools::conformedNcBoundaryField(bF)
    );

    // Replace every original patch field with a conformed patch field
    // containing the conformed orig and non-conformal fields
    forAll(origPatchIndices, i)
    {
        const label origPatchi = origPatchIndices[i];
        const fvPatch& origFvp = fvbm[origPatchi];

        autoPtr<conformedFvsPatchField<Type>> pF
        (
            new conformedFvsPatchField<Type>
            (
                origFvp,
                iF,
                bF.set(origPatchi, nullptr),
                autoPtr<calculatedFvsPatchField<Type>>
                (
                    new calculatedFvsPatchField<Type>(origFvp, iF)
                )
            )
        );

        pF() == origBf[origPatchi];
        pF->origFieldPtr_() == origBf[origPatchi];
        pF->ncFieldPtr_() == ncBf[origPatchi];

        bF.set(origPatchi, pF.ptr());
    }
}


template<class Type>
void Foam::conformedFvsPatchField<Type>::unconform
(
    typename SurfaceField<Type>::Boundary& bF
)
{
    const DimensionedField<Type, surfaceMesh>& iF = bF[0].internalField();

    const fvBoundaryMesh& fvbm = iF.mesh().boundary();

    const labelList origPatchIndices =
        nonConformalBoundary::New(iF.mesh()).allOrigPatchIndices();

    // Extract the conformed orig and non-conformal boundary fields from
    // the stored conformed patch fields
    PtrList<fvsPatchField<Type>> origPFs(fvbm.size());
    PtrList<fvsPatchField<Type>> ncPFs(fvbm.size());
    forAll(origPatchIndices, i)
    {
        const label origPatchi = origPatchIndices[i];

        conformedFvsPatchField<Type>& cpF =
            refCast<conformedFvsPatchField<Type>>(bF[origPatchi]);

        // If the mesh has topo-changed then maintained surface fields should
        // have been mapped or re-interpolated. So, copy the value from the
        // base field into the original field.
        if (iF.mesh().topoChanged())
        {
            cpF.origFieldPtr_() = bF[origPatchi];
        }

        origPFs.set(origPatchi, cpF.origFieldPtr_.ptr());
        ncPFs.set(origPatchi, cpF.ncFieldPtr_.ptr());
    }
    forAll(origPFs, patchi)
    {
        if (origPFs.set(patchi)) continue;

        origPFs.set(patchi, bF.set(patchi, nullptr));
        ncPFs.set
        (
            patchi,
            fvsPatchField<Type>::New
            (
                calculatedFvsPatchField<Type>::typeName,
                fvbm[patchi],
                iF
            )
        );
    }

    // If the mesh has topo-changed then just use the original parts and leave
    // the non-conformal parts unset
    if (iF.mesh().topoChanged())
    {
        bF.transfer(origPFs);
    }
    // If the mesh has not topo-changed, then combine the conformed boundary
    // fields to create the non-conformal boundary field
    else
    {
        typename SurfaceField<Type>::Boundary origAndNcBf
        (
            iF,
            fvMeshStitcherTools::unconformedBoundaryField
            (
                typename SurfaceField<Type>::Boundary(fvbm, iF, ncPFs),
                typename SurfaceField<Type>::Boundary(fvbm, iF, origPFs)
            )
        );

        bF.transfer(origAndNcBf);

        bF = fvMeshStitcherTools::synchronisedBoundaryField(bF);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::conformedFvsPatchField<Type>::conformedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(p, iF)
{
    NotImplemented;
}


template<class Type>
Foam::conformedFvsPatchField<Type>::conformedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchField<Type>(p, iF, dict),
    origFieldPtr_
    (
        fvsPatchField<Type>::New(p, iF, dict.subDict("origField")).ptr()
    ),
    ncFieldPtr_
    (
        new calculatedFvsPatchField<Type>(p, iF, dict.subDict("ncField"))
    )
{}


template<class Type>
Foam::conformedFvsPatchField<Type>::conformedFvsPatchField
(
    const conformedFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fieldMapper& mapper
)
:
    fvsPatchField<Type>(ptf, p, iF, mapper),
    origFieldPtr_
    (
        fvsPatchField<Type>::New(ptf.origFieldPtr_(), p, iF, mapper).ptr()
    ),
    ncFieldPtr_
    (
        new calculatedFvsPatchField<Type>
        (
            ptf.ncFieldPtr_(),
            p,
            iF,
            mapper
        )
    )
{}


template<class Type>
Foam::conformedFvsPatchField<Type>::conformedFvsPatchField
(
    const conformedFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(ptf, iF),
    origFieldPtr_(ptf.origFieldPtr_->clone(iF).ptr()),
    ncFieldPtr_(new calculatedFvsPatchField<Type>(ptf.ncFieldPtr_(), iF))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::conformedFvsPatchField<Type>::map
(
    const fvsPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    fvsPatchField<Type>::map(ptf, mapper);

    if (isA<conformedFvsPatchField<Type>>(ptf))
    {
        const conformedFvsPatchField<Type>& cptf =
            refCast<const conformedFvsPatchField<Type>>(ptf);

        origFieldPtr_->map(cptf.origFieldPtr_(), mapper);
        ncFieldPtr_->map(cptf.ncFieldPtr_(), mapper);
    }
    else
    {
        origFieldPtr_->reset(ptf);
        ncFieldPtr_() == origFieldPtr_();
    }
}


template<class Type>
void Foam::conformedFvsPatchField<Type>::reset(const fvsPatchField<Type>& ptf)
{
    fvsPatchField<Type>::reset(ptf);

    if (isA<conformedFvsPatchField<Type>>(ptf))
    {
        const conformedFvsPatchField<Type>& cptf =
            refCast<const conformedFvsPatchField<Type>>(ptf);

        origFieldPtr_->reset(cptf.origFieldPtr_());
        ncFieldPtr_->reset(cptf.ncFieldPtr_());
    }
    else
    {
        origFieldPtr_->reset(ptf);
        ncFieldPtr_() == origFieldPtr_();
    }
}


template<class Type>
void Foam::conformedFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
    writeEntry(os, "value", *this);

    writeKeyword(os, "origField") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    origFieldPtr_->write(os);
    os  << decrIndent << indent << token::END_BLOCK << nl;

    writeKeyword(os, "ncField") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    ncFieldPtr_->write(os);
    os  << decrIndent << indent << token::END_BLOCK << nl;
}


// ************************************************************************* //

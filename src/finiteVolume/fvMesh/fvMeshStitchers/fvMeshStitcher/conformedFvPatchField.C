/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "conformedFvPatchField.H"
#include "fvMeshStitcherTools.H"
#include "nonConformalBoundary.H"
#include "nonConformalFvPatch.H"
#include "nonConformalErrorFvPatch.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "forwardFieldMapper.H"
#include "forwardOrAssignFieldMapper.H"
#include "reverseInterpolativeFieldMapper.H"
#include "setSizeAndZeroFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * Private Static Member Functions * * * * * * * * * * * //

template<class Type>
Foam::labelList Foam::conformedFvPatchField<Type>::ncOrigNcField
(
    const fvBoundaryMesh& fvbm
)
{
    labelList result(fvbm.size(), -1);
    labelList origNcPatchCount(fvbm.size(), 0);

    // Global non-conformal patches
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;
        if (isA<processorCyclicFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        result[ncPatchi] = origNcPatchCount[ncFvp.origPatchIndex()] ++;
    }

    // Non-conformal processor cyclics
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;
        if (!isA<processorCyclicFvPatch>(fvp)) continue;

        const nonConformalProcessorCyclicFvPatch& ncpcFvp =
            refCast<const nonConformalProcessorCyclicFvPatch>(fvp);

        result[ncPatchi] = result[ncpcFvp.referPatchIndex()];
    }

    return result;
}


// * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * * //

template<class Type>
Foam::conformedFvPatchField<Type>::conformedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    autoPtr<fvPatchField<Type>>&& origFieldPtr,
    PtrList<fvPatchField<Type>>&& ncFieldPtrs,
    PtrList<scalarField>&& ncCoverages
)
:
    fvPatchField<Type>(p, iF),
    origFieldPtr_(origFieldPtr),
    ncFieldPtrs_(ncFieldPtrs),
    ncCoverages_(ncCoverages)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::labelList Foam::conformedFvPatchField<Type>::ncPatchIndices() const
{
    const fvPatch& origFvp = this->patch();
    const fvBoundaryMesh& fvbm = origFvp.boundaryMesh();

    // Determine the indices of the non-conformal patches associated with this
    // original patch. This assumes that all meshes that we might be mapping
    // between have their (non-processor) patches in the same order.
    DynamicList<label> result(ncFieldPtrs_.size());
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;
        if (isA<processorCyclicFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        if (ncFvp.origPatchName() != origFvp.name()) continue;

        result.append(ncPatchi);
    }

    // Check the number of identified indices corresponds to the number of
    // stored patch fields
    if (result.size() != ncFieldPtrs_.size())
    {
        const DimensionedField<Type, volMesh>& iF = this->internalField();

        FatalErrorInFunction
            << "Conformed field " << (isNull(iF) ? "" : iF.name() + ' ')
            << "on patch " << origFvp.name() << " has " << ncFieldPtrs_.size()
            << " non-conformal fields stored, but the patch has "
            << result.size() << " associated non-conformal patches"
            << exit(FatalError);
    }

    return result;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::conformedFvPatchField<Type>::conform
(
    typename VolField<Type>::Boundary& bF
)
{
    const DimensionedField<Type, volMesh>& iF = bF[0].internalField();

    const fvBoundaryMesh& fvbm = iF.mesh().boundary();

    const labelList origPatchIndices =
        nonConformalBoundary::New(iF.mesh()).allOrigPatchIndices();

    const labelList ncOrigNcField =
        conformedFvPatchField<Type>::ncOrigNcField(fvbm);

    // Replace every patch field on an orig boundary with a conformed patch
    // field, and move the previous orig boundary field into the member data of
    // the conformed patch field
    forAll(origPatchIndices, i)
    {
        const label origPatchi = origPatchIndices[i];
        const fvPatch& origFvp = fvbm[origPatchi];

        autoPtr<conformedFvPatchField<Type>> pF
        (
            new conformedFvPatchField<Type>
            (
                origFvp,
                iF,
                bF.set(origPatchi, nullptr),
                PtrList<fvPatchField<Type>>(),
                PtrList<scalarField>()
            )
        );

        pF() == pF->origFieldPtr_();

        bF.set(origPatchi, pF.ptr());
    }

    // Copy every (non-processor) non-conformal patch field into the conformed
    // patch field on the associated orig patch.
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;
        if (isA<processorCyclicFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        const label origPatchi = ncFvp.origPatchIndex();
        const fvPatch& origFvp = fvbm[origPatchi];

        conformedFvPatchField<Type>& cpF =
            refCast<conformedFvPatchField<Type>>(bF[origPatchi]);

        const label n =
            max(cpF.ncFieldPtrs_.size(), ncOrigNcField[ncPatchi] + 1);

        cpF.ncFieldPtrs_.resize(n);
        cpF.ncFieldPtrs_.set
        (
            ncOrigNcField[ncPatchi],
            fvPatchField<Type>::New
            (
                bF[ncPatchi],
                fvp,
                iF,
                setSizeAndZeroFieldMapper(origFvp.size())
            )
        );

        cpF.ncCoverages_.resize(n);
        cpF.ncCoverages_.set
        (
            ncOrigNcField[ncPatchi],
            scalarField(origFvp.size(), 0)
        );
    }

    // Compute the coverage for every non-conformal patch
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        const label origPatchi = ncFvp.origPatchIndex();
        const fvPatch& origFvp = fvbm[origPatchi];

        const labelList ncOrigPatchFace =
            ncFvp.polyFaces() - origFvp.start();

        const scalarField origNcMagSf
        (
            fvMeshStitcherTools::fieldRMapSum
            (
                ncFvp.patch().magSf(),
                origFvp.size(),
                ncOrigPatchFace
            )
        );

        conformedFvPatchField<Type>& cpF =
            refCast<conformedFvPatchField<Type>>(bF[origPatchi]);

        cpF.ncCoverages_[ncOrigNcField[ncPatchi]] +=
            origNcMagSf/origFvp.magSf();
    }

    // Map the data from the non-conformal patch fields to the conformed fields
    // by doing an area-weighted average
    forAll(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        const label origPatchi = ncFvp.origPatchIndex();
        const fvPatch& origFvp = fvbm[origPatchi];

        const labelList ncOrigPatchFace =
            ncFvp.polyFaces() - origFvp.start();

        conformedFvPatchField<Type>& cpF =
            refCast<conformedFvPatchField<Type>>(bF[origPatchi]);

        cpF.ncFieldPtrs_[ncOrigNcField[ncPatchi]].map
        (
            bF[ncPatchi],
            reverseInterpolativeFieldMapper
            (
                ncOrigPatchFace,
                ncFvp.patch().magSf()
               /max
                (
                    fvMeshStitcherTools::fieldMap
                    (
                        cpF.ncCoverages_[ncOrigNcField[ncPatchi]]
                       *origFvp.magSf(),
                        ncOrigPatchFace
                    ),
                    rootVSmall
                )
            )
        );
    }
}


template<class Type>
void Foam::conformedFvPatchField<Type>::unconform
(
    typename VolField<Type>::Boundary& bF
)
{
    const DimensionedField<Type, volMesh>& iF = bF[0].internalField();

    const fvBoundaryMesh& fvbm = iF.mesh().boundary();

    const labelList origPatchIndices =
        nonConformalBoundary::New(iF.mesh()).allOrigPatchIndices();

    const labelList ncOrigNcField =
        conformedFvPatchField<Type>::ncOrigNcField(fvbm);

    // Map out of the conformed fields into the non-conformal patches
    forAllReverse(fvbm, ncPatchi)
    {
        const fvPatch& fvp = fvbm[ncPatchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const nonConformalFvPatch& ncFvp =
            refCast<const nonConformalFvPatch>(fvp);

        const label origPatchi = ncFvp.origPatchIndex();
        const fvPatch& origFvp = fvbm[origPatchi];

        conformedFvPatchField<Type>& cpF =
            refCast<conformedFvPatchField<Type>>(bF[origPatchi]);

        labelList ncOrigPatchFace =
            ncFvp.polyFaces() - origFvp.start();

        forAll(ncOrigPatchFace, ncPatchFacei)
        {
            const label origPatchFacei = ncOrigPatchFace[ncPatchFacei];

            if (cpF.ncCoverages_[ncOrigNcField[ncPatchi]][origPatchFacei] == 0)
            {
                ncOrigPatchFace[ncPatchFacei] = -1;
            }
        }

        // If this is a calculated condition, then set unmapped face values to
        // the value in the adjacent cell. This is a tolerable approximation
        // and prevents errors in calculations that occur before the values are
        // re-computed as part of the solution.
        if (isType<calculatedFvPatchField<Type>>(bF[ncPatchi]))
        {
            bF[ncPatchi] = bF[ncPatchi].patchInternalField();

            bF[ncPatchi].map
            (
                cpF.ncFieldPtrs_[ncOrigNcField[ncPatchi]],
                forwardFieldMapper(ncOrigPatchFace)
            );
        }
        // Otherwise let the boundary condition determine what an appropriate
        // value is for an unmapped face
        else
        {
            bF[ncPatchi].map
            (
                cpF.ncFieldPtrs_[ncOrigNcField[ncPatchi]],
                forwardOrAssignPatchFieldMapper(bF[ncPatchi], ncOrigPatchFace)
            );
        }
    }

    // Move the orig patch fields into the boundary field
    forAll(origPatchIndices, i)
    {
        const label origPatchi = origPatchIndices[i];

        conformedFvPatchField<Type>& cpF =
            refCast<conformedFvPatchField<Type>>(bF[origPatchi]);

        bF.set(origPatchi, cpF.origFieldPtr_.ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::conformedFvPatchField<Type>::conformedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF)
{
    NotImplemented;
}


template<class Type>
Foam::conformedFvPatchField<Type>::conformedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    origFieldPtr_
    (
        fvPatchField<Type>::New(p, iF, dict.subDict("origField")).ptr()
    ),
    ncFieldPtrs_(dict.subDict("ncFields").size()),
    ncCoverages_(dict.subDict("ncFields").size())
{
    label i = 0;

    forAllConstIter(dictionary, dict.subDict("ncFields"), iter)
    {
        const word& ncPName = iter().keyword();
        const fvPatch& ncP = p.boundaryMesh()[ncPName];

        ncFieldPtrs_.set
        (
            i,
            fvPatchField<Type>::New(ncP, iF, iter().dict())
        );
        ncCoverages_.set
        (
            i,
            new scalarField("coverage", iter().dict(), p.size())
        );

        ++ i;
    }
}


template<class Type>
Foam::conformedFvPatchField<Type>::conformedFvPatchField
(
    const conformedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    origFieldPtr_
    (
        fvPatchField<Type>::New(ptf.origFieldPtr_(), p, iF, mapper).ptr()
    ),
    ncFieldPtrs_(ptf.ncFieldPtrs_.size()),
    ncCoverages_(ptf.ncCoverages_.size())
{
    const labelList ncPatchIndices(this->ncPatchIndices());

    forAll(ptf.ncFieldPtrs_, i)
    {
        ncFieldPtrs_.set
        (
            i,
            fvPatchField<Type>::New
            (
                ptf.ncFieldPtrs_[i],
                p.boundaryMesh()[ncPatchIndices[i]],
                iF,
                mapper
            )
        );
        ncCoverages_.set
        (
            i,
            mapper(ptf.ncCoverages_[i]).ptr()
        );
    }
}


template<class Type>
Foam::conformedFvPatchField<Type>::conformedFvPatchField
(
    const conformedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    origFieldPtr_(ptf.origFieldPtr_->clone(iF).ptr()),
    ncFieldPtrs_(ptf.ncFieldPtrs_, iF),
    ncCoverages_(ptf.ncCoverages_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::conformedFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    fvPatchField<Type>::map(ptf, mapper);

    if (isA<conformedFvPatchField<Type>>(ptf))
    {
        const conformedFvPatchField<Type>& cptf =
            refCast<const conformedFvPatchField<Type>>(ptf);

        origFieldPtr_->map(cptf.origFieldPtr_(), mapper);
        forAll(ncFieldPtrs_, i)
        {
            ncFieldPtrs_[i].map(cptf.ncFieldPtrs_[i], mapper);
            mapper(ncCoverages_[i], cptf.ncCoverages_[i]);
        }
    }
    else
    {
        NotImplemented;
    }
}


template<class Type>
void Foam::conformedFvPatchField<Type>::reset(const fvPatchField<Type>& ptf)
{
    fvPatchField<Type>::reset(ptf);

    if (isA<conformedFvPatchField<Type>>(ptf))
    {
        const conformedFvPatchField<Type>& cptf =
            refCast<const conformedFvPatchField<Type>>(ptf);

        origFieldPtr_->reset(cptf.origFieldPtr_());
        forAll(ncFieldPtrs_, i)
        {
            ncFieldPtrs_[i].reset(cptf.ncFieldPtrs_[i]);
            ncCoverages_[i] = cptf.ncCoverages_[i];
        }
    }
    else
    {
        NotImplemented;
    }
}


template<class Type>
void Foam::conformedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "value", *this);

    writeKeyword(os, "origField") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    origFieldPtr_->write(os);
    os  << decrIndent << indent << token::END_BLOCK << nl;

    writeKeyword(os, "ncFields") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    forAll(ncFieldPtrs_, i)
    {
        writeKeyword(os, ncFieldPtrs_[i].patch().name()) << nl;
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        ncFieldPtrs_[i].write(os);
        writeEntry(os, "coverage", ncCoverages_[i]);
        os  << decrIndent << indent << token::END_BLOCK << nl;
    }
    os  << decrIndent << indent << token::END_BLOCK << nl;
}


// ************************************************************************* //

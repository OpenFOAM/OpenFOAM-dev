/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "volFields.H"
#include "surfaceFields.H"
#include "calculatedFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "coupledFvPatchFields.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, facei)
    {
        intf[addr[facei]] += pf[facei];
    }
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, facei)
    {
        intf[addr[facei]] -= pf[facei];
    }
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(internalCoeffs_, patchi)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchi),
            internalCoeffs_[patchi].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    forAll(internalCoeffs_, patchi)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchi),
            cmptAv(internalCoeffs_[patchi]),
            diag
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::addBoundarySource
(
    Field<Type>& source,
    const bool couples
) const
{
    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];
        const Field<Type>& pbc = boundaryCoeffs_[patchi];

        if (!ptf.coupled())
        {
            addToInternalField(lduAddr().patchAddr(patchi), pbc, source);
        }
        else if (couples)
        {
            const tmp<Field<Type>> tpnf = ptf.patchNeighbourField();
            const Field<Type>& pnf = tpnf();

            const labelUList& addr = lduAddr().patchAddr(patchi);

            forAll(addr, facei)
            {
                source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
            }
        }
    }
}


template<class Type>
template<template<class> class ListType>
void Foam::fvMatrix<Type>::setValuesFromList
(
    const labelUList& cellLabels,
    const ListType<Type>& values
)
{
    const fvMesh& mesh = psi_.mesh();

    const cellList& cells = mesh.cells();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    scalarField& Diag = diag();
    Field<Type>& psi =
        const_cast
        <
            GeometricField<Type, fvPatchField, volMesh>&
        >(psi_).primitiveFieldRef();

    forAll(cellLabels, i)
    {
        const label celli = cellLabels[i];
        const Type& value = values[i];

        psi[celli] = value;
        source_[celli] = value*Diag[celli];

        if (symmetric() || asymmetric())
        {
            const cell& c = cells[celli];

            forAll(c, j)
            {
                const label facei = c[j];

                if (mesh.isInternalFace(facei))
                {
                    if (symmetric())
                    {
                        if (celli == own[facei])
                        {
                            source_[nei[facei]] -= upper()[facei]*value;
                        }
                        else
                        {
                            source_[own[facei]] -= upper()[facei]*value;
                        }

                        upper()[facei] = 0.0;
                    }
                    else
                    {
                        if (celli == own[facei])
                        {
                            source_[nei[facei]] -= lower()[facei]*value;
                        }
                        else
                        {
                            source_[own[facei]] -= upper()[facei]*value;
                        }

                        upper()[facei] = 0.0;
                        lower()[facei] = 0.0;
                    }
                }
                else
                {
                    label patchi = mesh.boundaryMesh().whichPatch(facei);

                    if (internalCoeffs_[patchi].size())
                    {
                        label patchFacei =
                            mesh.boundaryMesh()[patchi].whichFace(facei);

                        internalCoeffs_[patchi][patchFacei] =
                            Zero;

                        boundaryCoeffs_[patchi][patchFacei] =
                            Zero;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvMatrix<Type>::fvMatrix
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), Zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing fvMatrix<Type> for field " << psi_.name() << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchi)
    {
        internalCoeffs_.set
        (
            patchi,
            new Field<Type>
            (
                psi.mesh().boundary()[patchi].size(),
                Zero
            )
        );

        boundaryCoeffs_.set
        (
            patchi,
            new Field<Type>
            (
                psi.mesh().boundary()[patchi].size(),
                Zero
            )
        );
    }

    // Update the boundary coefficients of psi without changing its event No.
    GeometricField<Type, fvPatchField, volMesh>& psiRef =
       const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

    label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryFieldRef().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const fvMatrix<Type>& fvm)
:
    tmp<fvMatrix<Type>>::refCount(),
    lduMatrix(fvm),
    psi_(fvm.psi_),
    dimensions_(fvm.dimensions_),
    source_(fvm.source_),
    internalCoeffs_(fvm.internalCoeffs_),
    boundaryCoeffs_(fvm.boundaryCoeffs_),
    faceFluxCorrectionPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction
            << "Copying fvMatrix<Type> for field " << psi_.name() << endl;
    }

    if (fvm.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            *(fvm.faceFluxCorrectionPtr_)
        );
    }
}


template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const tmp<fvMatrix<Type>>& tfvm)
:
    lduMatrix
    (
        const_cast<fvMatrix<Type>&>(tfvm()),
        tfvm.isTmp()
    ),
    psi_(tfvm().psi_),
    dimensions_(tfvm().dimensions_),
    source_
    (
        const_cast<fvMatrix<Type>&>(tfvm()).source_,
        tfvm.isTmp()
    ),
    internalCoeffs_
    (
        const_cast<fvMatrix<Type>&>(tfvm()).internalCoeffs_,
        tfvm.isTmp()
    ),
    boundaryCoeffs_
    (
        const_cast<fvMatrix<Type>&>(tfvm()).boundaryCoeffs_,
        tfvm.isTmp()
    ),
    faceFluxCorrectionPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction
            << "Copying fvMatrix<Type> for field " << psi_.name() << endl;
    }

    if (tfvm().faceFluxCorrectionPtr_)
    {
        if (tfvm.isTmp())
        {
            faceFluxCorrectionPtr_ = tfvm().faceFluxCorrectionPtr_;
            tfvm().faceFluxCorrectionPtr_ = nullptr;
        }
        else
        {
            faceFluxCorrectionPtr_ = new
                GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    *(tfvm().faceFluxCorrectionPtr_)
                );
        }
    }

    tfvm.clear();
}


template<class Type>
Foam::fvMatrix<Type>::fvMatrix
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    Istream& is
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(is),
    source_(is),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing fvMatrix<Type> for field " << psi_.name() << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchi)
    {
        internalCoeffs_.set
        (
            patchi,
            new Field<Type>
            (
                psi.mesh().boundary()[patchi].size(),
                Zero
            )
        );

        boundaryCoeffs_.set
        (
            patchi,
            new Field<Type>
            (
                psi.mesh().boundary()[patchi].size(),
                Zero
            )
        );
    }

}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvMatrix<Type>::clone() const
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>(*this)
    );
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fvMatrix<Type>::~fvMatrix()
{
    if (debug)
    {
        InfoInFunction
            << "Destroying fvMatrix<Type> for field " << psi_.name() << endl;
    }

    if (faceFluxCorrectionPtr_)
    {
        delete faceFluxCorrectionPtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const UList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}


template<class Type>
void Foam::fvMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const UIndirectList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}


template<class Type>
void Foam::fvMatrix<Type>::setReference
(
    const label celli,
    const Type& value,
    const bool forceReference
)
{
    if ((forceReference || psi_.needReference()) && celli >= 0)
    {
        source()[celli] += diag()[celli]*value;
        diag()[celli] += diag()[celli];
    }
}


template<class Type>
void Foam::fvMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    if (debug)
    {
        InfoInFunction
            << "Relaxing " << psi_.name() << " by " << alpha << endl;
    }

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), 0.0);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchi);
            Field<Type>& iCoeffs = internalCoeffs_[patchi];

            if (ptf.coupled())
            {
                const Field<Type>& pCoeffs = boundaryCoeffs_[patchi];

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                forAll(pa, face)
                {
                    D[pa[face]] += component(iCoeffs[face], 0);
                    sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
                }
            }
            else
            {
                // For non-coupled boundaries add the maximum magnitude diagonal
                // contribution to ensure stability
                forAll(pa, face)
                {
                    D[pa[face]] += cmptMax(cmptMag(iCoeffs[face]));
                }
            }
        }
    }


    if (debug)
    {
        // Calculate amount of non-dominance.
        label nNon = 0;
        scalar maxNon = 0.0;
        scalar sumNon = 0.0;
        forAll(D, celli)
        {
            scalar d = (sumOff[celli] - D[celli])/mag(D[celli]);

            if (d > 0)
            {
                nNon++;
                maxNon = max(maxNon, d);
                sumNon += d;
            }
        }

        reduce(nNon, sumOp<label>(), UPstream::msgType(), psi_.mesh().comm());
        reduce
        (
            maxNon,
            maxOp<scalar>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );
        reduce
        (
            sumNon,
            sumOp<scalar>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );
        sumNon /= returnReduce
        (
            D.size(),
            sumOp<label>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );

        InfoInFunction
            << "Matrix dominance test for " << psi_.name() << nl
            << "    number of non-dominant cells   : " << nNon << nl
            << "    maximum relative non-dominance : " << maxNon << nl
            << "    average relative non-dominance : " << sumNon << nl
            << endl;
    }


    // Ensure the matrix is diagonally dominant...
    // Assumes that the central coefficient is positive and ensures it is
    forAll(D, celli)
    {
        D[celli] = max(mag(D[celli]), sumOff[celli]);
    }

    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchi);
            Field<Type>& iCoeffs = internalCoeffs_[patchi];

            if (ptf.coupled())
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= component(iCoeffs[face], 0);
                }
            }
            else
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= cmptMin(iCoeffs[face]);
                }
            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.primitiveField();
}


template<class Type>
void Foam::fvMatrix<Type>::relax()
{
    if
    (
        psi_.mesh().data::template lookupOrDefault<bool>
        ("finalIteration", false)
     && psi_.mesh().relaxEquation(psi_.name() + "Final")
    )
    {
        relax(psi_.mesh().equationRelaxationFactor(psi_.name() + "Final"));
    }
    else if (psi_.mesh().relaxEquation(psi_.name()))
    {
        relax(psi_.mesh().equationRelaxationFactor(psi_.name()));
    }
}


template<class Type>
void Foam::fvMatrix<Type>::boundaryManipulate
(
    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& bFields
)
{
    forAll(bFields, patchi)
    {
        bFields[patchi].manipulateMatrix(*this);
    }
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::fvMatrix<Type>::D() const
{
    tmp<scalarField> tdiag(new scalarField(diag()));
    addCmptAvBoundaryDiag(tdiag.ref());
    return tdiag;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvMatrix<Type>::DD() const
{
    tmp<Field<Type>> tdiag(pTraits<Type>::one*diag());

    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (!ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchAddr(patchi),
                internalCoeffs_[patchi],
                tdiag.ref()
            );
        }
    }

    return tdiag;
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Type>::A() const
{
    tmp<volScalarField> tAphi
    (
        volScalarField::New
        (
            "A("+psi_.name()+')',
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimVol,
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    tAphi.ref().primitiveFieldRef() = D()/psi_.mesh().V();
    tAphi.ref().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvMatrix<Type>::H() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tHphi
    (
        GeometricField<Type, fvPatchField, volMesh>::New
        (
            "H("+psi_.name()+')',
            psi_.mesh(),
            dimensions_/dimVol,
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Hphi = tHphi.ref();

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt(psi_.primitiveField().component(cmpt));

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        Hphi.primitiveFieldRef().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    Hphi.primitiveFieldRef() += lduMatrix::H(psi_.primitiveField()) + source_;
    addBoundarySource(Hphi.primitiveFieldRef());

    Hphi.primitiveFieldRef() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    typename Type::labelType validComponents
    (
        psi_.mesh().template validComponents<Type>()
    );

    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1)
        {
            Hphi.replace
            (
                cmpt,
                dimensionedScalar(Hphi.dimensions(), 0)
            );
        }
    }

    return tHphi;
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Type>::H1() const
{
    tmp<volScalarField> tH1
    (
        volScalarField::New
        (
            "H(1)",
            psi_.mesh(),
            dimensions_/(dimVol*psi_.dimensions()),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& H1_ = tH1.ref();

    H1_.primitiveFieldRef() = lduMatrix::H1();

    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchAddr(patchi),
                boundaryCoeffs_[patchi].component(0),
                H1_
            );
        }
    }

    H1_.primitiveFieldRef() /= psi_.mesh().V();
    H1_.correctBoundaryConditions();

    return tH1;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvMatrix<Type>::
flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorInFunction
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary"
               " of fvSchemes."
            << abort(FatalError);
    }

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfieldFlux
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            "flux("+psi_.name()+')',
            psi_.mesh(),
            dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& fieldFlux =
        tfieldFlux.ref();

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        fieldFlux.primitiveFieldRef().replace
        (
            cmpt,
            lduMatrix::faceH(psi_.primitiveField().component(cmpt))
        );
    }

    FieldField<Field, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchi)
    {
        InternalContrib[patchi] =
            cmptMultiply
            (
                InternalContrib[patchi],
                psi_.boundaryField()[patchi].patchInternalField()
            );
    }

    FieldField<Field, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchi)
    {
        if (psi_.boundaryField()[patchi].coupled())
        {
            NeighbourContrib[patchi] =
                cmptMultiply
                (
                    NeighbourContrib[patchi],
                    psi_.boundaryField()[patchi].patchNeighbourField()
                );
        }
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& ffbf = fieldFlux.boundaryFieldRef();

    forAll(ffbf, patchi)
    {
        ffbf[patchi] = InternalContrib[patchi] - NeighbourContrib[patchi];
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::operator=(const fvMatrix<Type>& fvmv)
{
    if (this == &fvmv)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(fvmv.psi_))
    {
        FatalErrorInFunction
            << "different fields"
            << abort(FatalError);
    }

    dimensions_ = fvmv.dimensions_;
    lduMatrix::operator=(fvmv);
    source_ = fvmv.source_;
    internalCoeffs_ = fvmv.internalCoeffs_;
    boundaryCoeffs_ = fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator=(const tmp<fvMatrix<Type>>& tfvmv)
{
    operator=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ += fvmv.dimensions_;
    lduMatrix::operator+=(fvmv);
    source_ += fvmv.source_;
    internalCoeffs_ += fvmv.internalCoeffs_;
    boundaryCoeffs_ += fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            *fvmv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=(const tmp<fvMatrix<Type>>& tfvmv)
{
    operator+=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "-=");

    dimensions_ -= fvmv.dimensions_;
    lduMatrix::operator-=(fvmv);
    source_ -= fvmv.source_;
    internalCoeffs_ -= fvmv.internalCoeffs_;
    boundaryCoeffs_ -= fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (-*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=(const tmp<fvMatrix<Type>>& tfvmv)
{
    operator-=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().V()*su.field();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().V()*su.field();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= psi().mesh().V()*su;
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += psi().mesh().V()*su;
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const zero&
)
{}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const zero&
)
{}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const volScalarField::Internal& dsf
)
{
    dimensions_ *= dsf.dimensions();
    lduMatrix::operator*=(dsf.field());
    source_ *= dsf.field();

    forAll(boundaryCoeffs_, patchi)
    {
        scalarField pisf
        (
            dsf.mesh().boundary()[patchi].patchInternalField(dsf.field())
        );

        internalCoeffs_[patchi] *= pisf;
        boundaryCoeffs_[patchi] *= pisf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorInFunction
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const tmp<volScalarField::Internal>& tdsf
)
{
    operator*=(tdsf());
    tdsf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const tmp<volScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator/=
(
    const volScalarField::Internal& dsf
)
{
    dimensions_ /= dsf.dimensions();
    lduMatrix::operator/=(dsf.field());
    source_ /= dsf.field();

    forAll(boundaryCoeffs_, patchi)
    {
        scalarField pisf
        (
            dsf.mesh().boundary()[patchi].patchInternalField(dsf.field())
        );

        internalCoeffs_[patchi] /= pisf;
        boundaryCoeffs_[patchi] /= pisf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorInFunction
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator/=
(
    const tmp<volScalarField::Internal>& tdsf
)
{
    operator/=(tdsf());
    tdsf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator/=
(
    const tmp<volScalarField>& tvsf
)
{
    operator/=(tvsf());
    tvsf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator/=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ /= ds.dimensions();
    lduMatrix::operator/=(ds.value());
    source_ /= ds.value();
    internalCoeffs_ /= ds.value();
    boundaryCoeffs_ /= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ /= ds.value();
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm1,
    const fvMatrix<Type>& fvm2,
    const char* op
)
{
    if (&fvm1.psi() != &fvm2.psi())
    {
        FatalErrorInFunction
            << "incompatible fields for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << "] "
            << op
            << " [" << fvm2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fvm1.dimensions() != fvm2.dimensions())
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << fvm1.dimensions()/dimVolume << " ] "
            << op
            << " [" << fvm2.psi().name() << fvm2.dimensions()/dimVolume << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm,
    const DimensionedField<Type, volMesh>& df,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != df.dimensions())
    {
        FatalErrorInFunction
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << df.name() << df.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != dt.dimensions())
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::correction
(
    const fvMatrix<Type>& A
)
{
    tmp<Foam::fvMatrix<Type>> tAcorr = A - (A & A.psi());

    // Delete the faceFluxCorrection from the correction matrix
    // as it does not have a clear meaning or purpose
    deleteDemandDrivenData(tAcorr.ref().faceFluxCorrectionPtr());

    return tAcorr;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::correction
(
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<Foam::fvMatrix<Type>> tAcorr = tA - (tA() & tA().psi());

    // Delete the faceFluxCorrection from the correction matrix
    // as it does not have a clear meaning or purpose
    deleteDemandDrivenData(tAcorr.ref().faceFluxCorrectionPtr());

    return tAcorr;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += A.psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tC().psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const zero&
)
{
    return A;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const zero&
)
{
    return tA;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<fvMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<fvMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const volScalarField::Internal& dsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField::Internal>& tdsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const volScalarField::Internal& dsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField::Internal>& tdsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= ds;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const fvMatrix<Type>& A,
    const volScalarField::Internal& dsf
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() /= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const fvMatrix<Type>& A,
    const tmp<volScalarField::Internal>& tdsf
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() /= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const fvMatrix<Type>& A,
    const tmp<volScalarField>& tvsf
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() /= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const tmp<fvMatrix<Type>>& tA,
    const volScalarField::Internal& dsf
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() /= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<volScalarField::Internal>& tdsf
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() /= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<volScalarField>& tvsf
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() /= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const fvMatrix<Type>& A,
    const dimensioned<scalar>& ds
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() /= ds;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator/
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<scalar>& ds
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() /= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMphi
    (
        GeometricField<Type, fvPatchField, volMesh>::New
        (
            "M&" + psi.name(),
            psi.mesh(),
            M.dimensions()/dimVol,
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Mphi = tMphi.ref();

    // Loop over field components
    if (M.hasDiag())
    {
        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            scalarField psiCmpt(psi.field().component(cmpt));
            scalarField boundaryDiagCmpt(M.diag());
            M.addBoundaryDiag(boundaryDiagCmpt, cmpt);
            Mphi.primitiveFieldRef().replace(cmpt, -boundaryDiagCmpt*psiCmpt);
        }
    }
    else
    {
        Mphi.primitiveFieldRef() = Zero;
    }

    Mphi.primitiveFieldRef() += M.lduMatrix::H(psi.field()) + M.source();
    M.addBoundarySource(Mphi.primitiveFieldRef());

    Mphi.primitiveFieldRef() /= -psi.mesh().V();
    Mphi.correctBoundaryConditions();

    return tMphi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<DimensionedField<Type, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & psi;
    tM.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const tmp<DimensionedField<Type, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvMatrix<Type>& fvm)
{
    os  << static_cast<const lduMatrix&>(fvm) << nl
        << fvm.dimensions_ << nl
        << fvm.source_ << nl
        << fvm.internalCoeffs_ << nl
        << fvm.boundaryCoeffs_ << endl;

    os.check("Ostream& operator<<(Ostream&, fvMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "fvMatrixSolve.C"

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "zeroGradientFvPatchFields.H"
#include "coupledFvPatchFields.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
        FatalErrorIn
        (
            "fvMatrix<Type>::addToInternalField(const labelUList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] += pf[faceI];
    }
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2> >& tpf,
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
        FatalErrorIn
        (
            "fvMatrix<Type>::addToInternalField(const labelUList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] -= pf[faceI];
    }
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2> >& tpf,
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
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            internalCoeffs_[patchI].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            cmptAv(internalCoeffs_[patchI]),
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
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];
        const Field<Type>& pbc = boundaryCoeffs_[patchI];

        if (!ptf.coupled())
        {
            addToInternalField(lduAddr().patchAddr(patchI), pbc, source);
        }
        else if (couples)
        {
            tmp<Field<Type> > tpnf = ptf.patchNeighbourField();
            const Field<Type>& pnf = tpnf();

            const labelUList& addr = lduAddr().patchAddr(patchI);

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
        >(psi_).internalField();

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
                            pTraits<Type>::zero;

                        boundaryCoeffs_[patchi][patchFacei] =
                            pTraits<Type>::zero;
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
    source_(psi.size(), pTraits<Type>::zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>(GeometricField<Type, fvPatchField, volMesh>&,"
               " const dimensionSet&) : "
               "constructing fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

    // Update the boundary coefficients of psi without changing its event No.
    GeometricField<Type, fvPatchField, volMesh>& psiRef =
       const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

    label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryField().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const fvMatrix<Type>& fvm)
:
    refCount(),
    lduMatrix(fvm),
    psi_(fvm.psi_),
    dimensions_(fvm.dimensions_),
    source_(fvm.source_),
    internalCoeffs_(fvm.internalCoeffs_),
    boundaryCoeffs_(fvm.boundaryCoeffs_),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::fvMatrix(const fvMatrix<Type>&) : "
            << "copying fvMatrix<Type> for field " << psi_.name()
            << endl;
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


#ifdef ConstructFromTmp
template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const tmp<fvMatrix<Type> >& tfvm)
:
    refCount(),
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
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::fvMatrix(const tmp<fvMatrix<Type> >&) : "
            << "copying fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (tfvm().faceFluxCorrectionPtr_)
    {
        if (tfvm.isTmp())
        {
            faceFluxCorrectionPtr_ = tfvm().faceFluxCorrectionPtr_;
            tfvm().faceFluxCorrectionPtr_ = NULL;
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
#endif


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
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>"
               "(GeometricField<Type, fvPatchField, volMesh>&, Istream&) : "
               "constructing fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

}


template<class Type>
Foam::fvMatrix<Type>::~fvMatrix()
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::~fvMatrix<Type>() : "
            << "destroying fvMatrix<Type> for field " << psi_.name()
            << endl;
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
        InfoIn("fvMatrix<Type>::relax(const scalar alpha)")
            << "Relaxing " << psi_.name() << " by " << alpha
            << endl;
    }

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), 0.0);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                const Field<Type>& pCoeffs = boundaryCoeffs_[patchI];

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

        InfoIn("fvMatrix<Type>::relax(const scalar alpha)")
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
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

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
    S += (D - D0)*psi_.internalField();
}


template<class Type>
void Foam::fvMatrix<Type>::relax()
{
    word name = psi_.select
    (
        psi_.mesh().data::template lookupOrDefault<bool>
        ("finalIteration", false)
    );

    if (psi_.mesh().relaxEquation(name))
    {
        relax(psi_.mesh().equationRelaxationFactor(name));
    }
}


template<class Type>
void Foam::fvMatrix<Type>::boundaryManipulate
(
    typename GeometricField<Type, fvPatchField, volMesh>::
        GeometricBoundaryField& bFields
)
{
    forAll(bFields, patchI)
    {
        bFields[patchI].manipulateMatrix(*this);
    }
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::fvMatrix<Type>::D() const
{
    tmp<scalarField> tdiag(new scalarField(diag()));
    addCmptAvBoundaryDiag(tdiag());
    return tdiag;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fvMatrix<Type>::DD() const
{
    tmp<Field<Type> > tdiag(pTraits<Type>::one*diag());

    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (!ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchAddr(patchI),
                internalCoeffs_[patchI],
                tdiag()
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
        new volScalarField
        (
            IOobject
            (
                "A("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tAphi().internalField() = D()/psi_.mesh().V();
    tAphi().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::fvMatrix<Type>::H() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tHphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Hphi = tHphi();

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt(psi_.internalField().component(cmpt));

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        Hphi.internalField().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    Hphi.internalField() += lduMatrix::H(psi_.internalField()) + source_;
    addBoundarySource(Hphi.internalField());

    Hphi.internalField() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    typename Type::labelType validComponents
    (
        pow
        (
            psi_.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1)
        {
            Hphi.replace
            (
                cmpt,
                dimensionedScalar("0", Hphi.dimensions(), 0.0)
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
        new volScalarField
        (
            IOobject
            (
                "H(1)",
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/(dimVol*psi_.dimensions()),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& H1_ = tH1();

    H1_.internalField() = lduMatrix::H1();

    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchAddr(patchI),
                boundaryCoeffs_[patchI].component(0),
                H1_
            );
        }
    }

    H1_.internalField() /= psi_.mesh().V();
    H1_.correctBoundaryConditions();

    return tH1;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::fvMatrix<Type>::
flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorIn("fvMatrix<Type>::flux()")
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary"
               " of fvSchemes."
            << abort(FatalError);
    }

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfieldFlux
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& fieldFlux = tfieldFlux();

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        fieldFlux.internalField().replace
        (
            cmpt,
            lduMatrix::faceH(psi_.internalField().component(cmpt))
        );
    }

    FieldField<Field, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchI)
    {
        InternalContrib[patchI] =
            cmptMultiply
            (
                InternalContrib[patchI],
                psi_.boundaryField()[patchI].patchInternalField()
            );
    }

    FieldField<Field, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchI)
    {
        if (psi_.boundaryField()[patchI].coupled())
        {
            NeighbourContrib[patchI] =
                cmptMultiply
                (
                    NeighbourContrib[patchI],
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
        }
    }

    forAll(fieldFlux.boundaryField(), patchI)
    {
        fieldFlux.boundaryField()[patchI] =
            InternalContrib[patchI] - NeighbourContrib[patchI];
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
        FatalErrorIn("fvMatrix<Type>::operator=(const fvMatrix<Type>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(fvmv.psi_))
    {
        FatalErrorIn("fvMatrix<Type>::operator=(const fvMatrix<Type>&)")
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
void Foam::fvMatrix<Type>::operator=(const tmp<fvMatrix<Type> >& tfvmv)
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
void Foam::fvMatrix<Type>::operator+=(const tmp<fvMatrix<Type> >& tfvmv)
{
    operator+=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

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
void Foam::fvMatrix<Type>::operator-=(const tmp<fvMatrix<Type> >& tfvmv)
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
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
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
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
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
    const DimensionedField<scalar, volMesh>& dsf
)
{
    dimensions_ *= dsf.dimensions();
    lduMatrix::operator*=(dsf.field());
    source_ *= dsf.field();

    forAll(boundaryCoeffs_, patchI)
    {
        scalarField psf
        (
            dsf.mesh().boundary()[patchI].patchInternalField(dsf.field())
        );

        internalCoeffs_[patchI] *= psf;
        boundaryCoeffs_[patchI] *= psf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::operator*="
            "(const DimensionedField<scalar, volMesh>&)"
        )   << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const tmp<DimensionedField<scalar, volMesh> >& tdsf
)
{
    operator*=(tdsf());
    tdsf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const volScalarField& vsf
)
{
    dimensions_ *= vsf.dimensions();
    lduMatrix::operator*=(vsf.field());
    source_ *= vsf.field();

    forAll(vsf.boundaryField(), patchI)
    {
        const fvPatchScalarField& psf = vsf.boundaryField()[patchI];

        if (psf.coupled())
        {
            internalCoeffs_[patchI] *= psf.patchInternalField();
            boundaryCoeffs_[patchI] *= psf.patchNeighbourField();
        }
        else
        {
            internalCoeffs_[patchI] *= psf.patchInternalField();
            boundaryCoeffs_[patchI] *= psf;
        }
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::operator*="
            "(const DimensionedField<scalar, volMesh>&)"
        )   << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
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
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const fvMatrix<Type>&)"
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << "] "
            << op
            << " [" << fvm2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fvm1.dimensions() != fvm2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const fvMatrix<Type>&)"
        )   << "incompatible dimensions for operation "
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
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const GeometricField<Type, "
            "fvPatchField, volMesh>&)"
        )   <<  "incompatible dimensions for operation "
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
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const dimensioned<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::solverPerformance Foam::solve
(
    fvMatrix<Type>& fvm,
    const dictionary& solverControls
)
{
    return fvm.solve(solverControls);
}

template<class Type>
Foam::solverPerformance Foam::solve
(
    const tmp<fvMatrix<Type> >& tfvm,
    const dictionary& solverControls
)
{
    solverPerformance solverPerf =
        const_cast<fvMatrix<Type>&>(tfvm()).solve(solverControls);

    tfvm.clear();

    return solverPerf;
}


template<class Type>
Foam::solverPerformance Foam::solve(fvMatrix<Type>& fvm)
{
    return fvm.solve();
}

template<class Type>
Foam::solverPerformance Foam::solve(const tmp<fvMatrix<Type> >& tfvm)
{
    solverPerformance solverPerf =
        const_cast<fvMatrix<Type>&>(tfvm()).solve();

    tfvm.clear();

    return solverPerf;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::correction
(
    const fvMatrix<Type>& A
)
{
    tmp<Foam::fvMatrix<Type> > tAcorr = A - (A & A.psi());

    if
    (
        (A.hasUpper() || A.hasLower())
     && A.psi().mesh().fluxRequired(A.psi().name())
    )
    {
        tAcorr().faceFluxCorrectionPtr() = (-A.flux()).ptr();
    }

    return tAcorr;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::correction
(
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<Foam::fvMatrix<Type> > tAcorr = tA - (tA() & tA().psi());

    // Note the matrix coefficients are still that of matrix A
    const fvMatrix<Type>& A = tAcorr();

    if
    (
        (A.hasUpper() || A.hasLower())
     && A.psi().mesh().fluxRequired(A.psi().name())
    )
    {
        tAcorr().faceFluxCorrectionPtr() = (-A.flux()).ptr();
    }

    return tAcorr;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += A.psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tC().psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const zero&
)
{
    return A;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const zero&
)
{
    return tA;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<fvMatrix<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<fvMatrix<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const DimensionedField<scalar, volMesh>& dsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp< DimensionedField<scalar, volMesh> >& tdsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const DimensionedField<scalar, volMesh>& dsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp<DimensionedField<scalar, volMesh> >& tdsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const dimensioned<scalar>& ds,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= ds;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const dimensioned<scalar>& ds,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const fvMatrix<Type>& M,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "M&" + psi.name(),
                psi.instance(),
                psi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi.mesh(),
            M.dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Mphi = tMphi();

    // Loop over field components
    if (M.hasDiag())
    {
        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            scalarField psiCmpt(psi.field().component(cmpt));
            scalarField boundaryDiagCmpt(M.diag());
            M.addBoundaryDiag(boundaryDiagCmpt, cmpt);
            Mphi.internalField().replace(cmpt, -boundaryDiagCmpt*psiCmpt);
        }
    }
    else
    {
        Mphi.internalField() = pTraits<Type>::zero;
    }

    Mphi.internalField() += M.lduMatrix::H(psi.field()) + M.source();
    M.addBoundarySource(Mphi.internalField());

    Mphi.internalField() /= -psi.mesh().V();
    Mphi.correctBoundaryConditions();

    return tMphi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<DimensionedField<Type, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const tmp<fvMatrix<Type> >& tM,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & psi;
    tM.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const tmp<fvMatrix<Type> >& tM,
    const tmp<DimensionedField<Type, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const tmp<fvMatrix<Type> >& tM,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & tpsi();
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

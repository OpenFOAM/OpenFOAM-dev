/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "cyclicACMIFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{
    if (!isA<cyclicACMIFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField"
            "("
                "const cyclicACMIFvPatchField<Type>& ,"
                "const fvPatch&, "
                "const DimensionedField<Type, volMesh>&, "
                "const fvPatchFieldMapper&"
            ")"
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{
    if (!isA<cyclicACMIFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField"
            "("
                "const fvPatch&, "
                "const DimensionedField<Type, volMesh>&, "
                "const dictionary&"
            ")",
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::blocking);
    }
}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_)
{}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicACMIFvPatchField<Type>::coupled() const
{
    return cyclicACMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::cyclicACMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();
    const labelUList& nbrFaceCellsCoupled =
        cyclicACMIPatch_.cyclicACMIPatch().neighbPatch().faceCells();
    const labelUList& faceCellsNonOverlap =
        cyclicACMIPatch_.cyclicACMIPatch().nonOverlapPatch().faceCells();

    Field<Type> pnfCoupled(iField, nbrFaceCellsCoupled);
    Field<Type> pfNonOverlap(iField, faceCellsNonOverlap);

    tmp<Field<Type> > tpnf
    (
        new Field<Type>
        (
            cyclicACMIPatch_.interpolate
            (
                pnfCoupled,
                pfNonOverlap
            )
        )
    );

    if (doTransform())
    {
        tpnf() = transform(forwardT(), tpnf());
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicACMIFvPatchField<Type>&
Foam::cyclicACMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->internalField()
        );

    return refCast<const cyclicACMIFvPatchField<Type> >
    (
        fld.boundaryField()[cyclicACMIPatch_.neighbPatchID()]
    );
}


template<class Type>
const Foam::fvPatchField<Type>&
Foam::cyclicACMIFvPatchField<Type>::nonOverlapPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->internalField()
        );

    return fld.boundaryField()[cyclicACMIPatch_.nonOverlapPatchID()];
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // note: only applying coupled contribution

    const labelUList& nbrFaceCellsCoupled =
        cyclicACMIPatch_.cyclicACMIPatch().neighbPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCellsCoupled);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    const labelUList& faceCells = cyclicACMIPatch_.faceCells();

    pnf = cyclicACMIPatch_.interpolate(pnf);

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    // note: only applying coupled contribution

    const labelUList& nbrFaceCellsCoupled =
        cyclicACMIPatch_.cyclicACMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCellsCoupled);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    const labelUList& faceCells = cyclicACMIPatch_.faceCells();

    pnf = cyclicACMIPatch_.interpolate(pnf);

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::cyclicACMIFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    // note: only applying coupled contribution
    return coupledFvPatchField<Type>::snGrad(deltaCoeffs);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateCoeffs()
{
    // update non-overlap patch - some will implement updateCoeffs, and
    // others will implement evaluate

    // scale neighbour field by (1 - mask)

    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();
    const fvPatchField<Type>& npf = nonOverlapPatchField();
    const_cast<fvPatchField<Type>&>(npf).updateCoeffs(1.0 - mask);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes comms
)
{
    // update non-overlap patch (if not already updated by updateCoeffs)

    // scale neighbour field by (1 - mask)

    fvPatchField<Type>& npf =
        const_cast<fvPatchField<Type>&>(nonOverlapPatchField());

    if (!npf.updated())
    {
        const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();

        npf.evaluate(comms);

        npf *= 1.0 - mask;
    }
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes comms
)
{
    // blend contributions from the coupled and non-overlap patches

    // neighbour patch field is updated via updateCoeffs or initEvaluate
    // and is already scaled by (1 - mask)
    const fvPatchField<Type>& npf = nonOverlapPatchField();

    coupledFvPatchField<Type>::evaluate(comms);
    const Field<Type>& cpf = *this;

    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();
    Field<Type>::operator=(mask*cpf + npf);

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::cyclicACMIFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    // note: do not blend based on mask field
    // - when applied this is scaled by the areas which are already scaled
    return coupledFvPatchField<Type>::valueInternalCoeffs(w);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::cyclicACMIFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    // note: do not blend based on mask field
    // - when applied this is scaled by the areas which are already scaled
    return coupledFvPatchField<Type>::valueBoundaryCoeffs(w);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::cyclicACMIFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    // note: do not blend based on mask field
    // - when applied this is scaled by the areas which are already scaled
    return coupledFvPatchField<Type>::gradientInternalCoeffs(deltaCoeffs);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::cyclicACMIFvPatchField<Type>::gradientInternalCoeffs() const
{
    // note: do not blend based on mask field
    // - when applied this is scaled by the areas which are already scaled
    return coupledFvPatchField<Type>::gradientInternalCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::cyclicACMIFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    // note: do not blend based on mask field
    // - when applied this is scaled by the areas which are already scaled
    return coupledFvPatchField<Type>::gradientBoundaryCoeffs(deltaCoeffs);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::cyclicACMIFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    // note: do not blend based on mask field
    // - when applied this is scaled by the areas which are already scaled
    return coupledFvPatchField<Type>::gradientBoundaryCoeffs();
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();

    // nothing to be done by the AMI, but re-direct to non-overlap patch
    // with non-overlap patch weights
    const fvPatchField<Type>& npf = nonOverlapPatchField();

    const_cast<fvPatchField<Type>&>(npf).manipulateMatrix(matrix, 1.0 - mask);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //

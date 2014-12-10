/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "cyclicPointPatchField.H"
#include "Swap.H"
#include "transformField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicPointPatch>(p))
{}


template<class Type>
Foam::cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicPointPatch>(p))
{
    if (!isType<cyclicPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicPointPatchField<Type>::cyclicPointPatchField\n"
            "(\n"
            "    const pointPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not cyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const cyclicPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicPointPatch>(p))
{
    if (!isType<cyclicPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicPointPatchField<Type>::cyclicPointPatchField\n"
            "(\n"
            "    const cyclicPointPatchField<Type>& ptf,\n"
            "    const pointPatch& p,\n"
            "    const DimensionedField<Type, pointMesh>& iF,\n"
            "    const pointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const cyclicPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::cyclicPointPatchField<Type>::swapAddSeparated
(
    const Pstream::commsTypes,
    Field<Type>& pField
) const
{
    // Get neighbouring pointPatch
    const cyclicPointPatch& nbrPatch = cyclicPatch_.neighbPatch();

    if (cyclicPatch_.cyclicPatch().owner())
    {
        // We inplace modify pField. To prevent the other side (which gets
        // evaluated at a later date) using already changed values we do
        // all swaps on the side that gets evaluated first.

        // Get neighbouring pointPatchField
        const GeometricField<Type, pointPatchField, pointMesh>& fld =
            refCast<const GeometricField<Type, pointPatchField, pointMesh> >
            (
                this->dimensionedInternalField()
            );

        const cyclicPointPatchField<Type>& nbr =
            refCast<const cyclicPointPatchField<Type> >
            (
                fld.boundaryField()[nbrPatch.index()]
            );


        Field<Type> pf(this->patchInternalField(pField));
        Field<Type> nbrPf(nbr.patchInternalField(pField));

        const edgeList& pairs = cyclicPatch_.transformPairs();

        if (doTransform())
        {
            // Transform both sides.
            forAll(pairs, pairi)
            {
                label pointi = pairs[pairi][0];
                label nbrPointi = pairs[pairi][1];

                Type tmp = pf[pointi];
                pf[pointi] = transform(forwardT()[0], nbrPf[nbrPointi]);
                nbrPf[nbrPointi] = transform(reverseT()[0], tmp);
            }
        }
        else
        {
            forAll(pairs, pairi)
            {
                Swap(pf[pairs[pairi][0]], nbrPf[pairs[pairi][1]]);
            }
        }
        this->addToInternalField(pField, pf);
        nbr.addToInternalField(pField, nbrPf);
    }
}


// ************************************************************************* //

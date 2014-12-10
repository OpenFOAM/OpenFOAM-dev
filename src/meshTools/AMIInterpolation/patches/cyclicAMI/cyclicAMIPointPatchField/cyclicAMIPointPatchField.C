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

#include "cyclicAMIPointPatchField.H"
#include "Swap.H"
#include "transformField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIPointPatch>(p)),
    ppiPtr_(NULL),
    nbrPpiPtr_(NULL)
{}


template<class Type>
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    cyclicAMIPatch_(refCast<const cyclicAMIPointPatch>(p)),
    ppiPtr_(NULL),
    nbrPpiPtr_(NULL)
{
    if (!isType<cyclicAMIPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField\n"
            "(\n"
            "    const pointPatch&,\n"
            "    const DimensionedField<Type, pointMesh>&,\n"
            "    const dictionary&\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not cyclicAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const cyclicAMIPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIPointPatch>(p)),
    ppiPtr_(NULL),
    nbrPpiPtr_(NULL)
{
    if (!isType<cyclicAMIPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField\n"
            "(\n"
            "    const cyclicAMIPointPatchField<Type>&,\n"
            "    const pointPatch&,\n"
            "    const DimensionedField<Type, pointMesh>&,\n"
            "    const pointPatchFieldMapper&\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const cyclicAMIPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    ppiPtr_(NULL),
    nbrPpiPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::cyclicAMIPointPatchField<Type>::swapAddSeparated
(
    const Pstream::commsTypes,
    Field<Type>& pField
) const
{
    if (cyclicAMIPatch_.cyclicAMIPatch().owner())
    {
        // We inplace modify pField. To prevent the other side (which gets
        // evaluated at a later date) using already changed values we do
        // all swaps on the side that gets evaluated first.

        // Get neighbouring pointPatch
        const cyclicAMIPointPatch& nbrPatch = cyclicAMIPatch_.neighbPatch();

        // Get neighbouring pointPatchField
        const GeometricField<Type, pointPatchField, pointMesh>& fld =
            refCast<const GeometricField<Type, pointPatchField, pointMesh> >
            (
                this->dimensionedInternalField()
            );

        const cyclicAMIPointPatchField<Type>& nbr =
            refCast<const cyclicAMIPointPatchField<Type> >
            (
                fld.boundaryField()[nbrPatch.index()]
            );


        Field<Type> ptFld(this->patchInternalField(pField));
        Field<Type> nbrPtFld(nbr.patchInternalField(pField));


        if (doTransform())
        {
            const tensor& forwardT = this->forwardT()[0];
            const tensor& reverseT = this->reverseT()[0];

            transform(ptFld, reverseT, ptFld);
            transform(nbrPtFld, forwardT, nbrPtFld);
        }

        // convert point field to face field, AMI interpolate, then
        // face back to point
        {
            // add neighbour side contribution to owner
            Field<Type> nbrFcFld(nbrPpi().pointToFaceInterpolate(nbrPtFld));

            // interpolate to owner
            if (cyclicAMIPatch_.cyclicAMIPatch().applyLowWeightCorrection())
            {
                Field<Type> fcFld(ppi().pointToFaceInterpolate(ptFld));

                nbrFcFld =
                    cyclicAMIPatch_.cyclicAMIPatch().interpolate
                    (
                        nbrFcFld,
                        fcFld
                    );
            }
            else
            {
                nbrFcFld =
                    cyclicAMIPatch_.cyclicAMIPatch().interpolate(nbrFcFld);
            }

            // add to internal field
            this->addToInternalField
            (
                pField,
                ppi().faceToPointInterpolate(nbrFcFld)()
            );
        }

        {
            // add owner side contribution to neighbour
            Field<Type> fcFld(ppi().pointToFaceInterpolate(ptFld));

            // interpolate to neighbour
            if (cyclicAMIPatch_.cyclicAMIPatch().applyLowWeightCorrection())
            {
                Field<Type> nbrFcFld(nbrPpi().pointToFaceInterpolate(nbrPtFld));

                fcFld =
                    cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().interpolate
                    (
                        fcFld,
                        nbrFcFld
                    );
            }
            else
            {
                fcFld =
                    cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().interpolate
                    (
                        fcFld
                    );
            }

            // add to internal field
            nbr.addToInternalField
            (
                pField,
                nbrPpi().faceToPointInterpolate(fcFld)()
            );
        }
    }
}


// ************************************************************************* //

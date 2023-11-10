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

#include "processorFvsPatchField.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p))
{}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvsPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p))
{}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isA<processorFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not processor type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const processorFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isA<processorFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const processorFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvsPatchField<Type>::~processorFvsPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorFvsPatchField<Type>::initPatchNeighbourField
(
    const Pstream::commsTypes commsType
) const
{
    if (Pstream::parRun())
    {
        if (commsType == Pstream::commsTypes::nonBlocking)
        {
            receiveBuf_.setSize(this->size());

            IPstream::read
            (
                commsType,
                procPatch_.neighbProcNo(),
                reinterpret_cast<char*>(receiveBuf_.begin()),
                receiveBuf_.byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }

        OPstream::write
        (
            commsType,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(this->begin()),
            this->byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::processorFvsPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes commsType
) const
{
    if (Pstream::parRun())
    {
        if (commsType != Pstream::commsTypes::nonBlocking)
        {
            receiveBuf_.setSize(this->size());

            IPstream::read
            (
                commsType,
                procPatch_.neighbProcNo(),
                reinterpret_cast<char*>(receiveBuf_.begin()),
                receiveBuf_.byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }

        procPatch_.transform().transform(receiveBuf_, receiveBuf_);

        return receiveBuf_;
    }
    else
    {
        return tmp<Field<Type>>(new Field<Type>());
    }
}


// ************************************************************************* //

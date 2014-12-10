/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "processorCyclicFvPatchField.H"
#include "processorCyclicFvPatch.H"
#include "demandDrivenData.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorCyclicFvPatchField<Type>::processorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    processorFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorCyclicFvPatch>(p))
{}


template<class Type>
Foam::processorCyclicFvPatchField<Type>::processorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    //coupledFvPatchField<Type>(p, iF, f),
    processorFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorCyclicFvPatch>(p))
{}


// Construct by mapping given processorCyclicFvPatchField<Type>
template<class Type>
Foam::processorCyclicFvPatchField<Type>::processorCyclicFvPatchField
(
    const processorCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    //coupledFvPatchField<Type>(ptf, p, iF, mapper),
    processorFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorCyclicFvPatch>(p))
{
    if (!isType<processorCyclicFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "processorCyclicFvPatchField<Type>::processorCyclicFvPatchField\n"
            "(\n"
            "    const processorCyclicFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorCyclicFvPatchField<Type>::processorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    //coupledFvPatchField<Type>(p, iF, dict),
    processorFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorCyclicFvPatch>(p))
{
    if (!isType<processorCyclicFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "processorCyclicFvPatchField<Type>::processorCyclicFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        WarningIn
        (
            "processorCyclicFvPatchField<Type>::processorCyclicFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "Scheduled communication with split cyclics not supported."
            << endl;
    }
}


template<class Type>
Foam::processorCyclicFvPatchField<Type>::processorCyclicFvPatchField
(
    const processorCyclicFvPatchField<Type>& ptf
)
:
    //processorLduInterfaceField(),
    //coupledFvPatchField<Type>(ptf),
    processorFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorCyclicFvPatch>(ptf.patch()))
{}


template<class Type>
Foam::processorCyclicFvPatchField<Type>::processorCyclicFvPatchField
(
    const processorCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    //coupledFvPatchField<Type>(ptf, iF),
    processorFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorCyclicFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorCyclicFvPatchField<Type>::~processorCyclicFvPatchField()
{}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "CloudAverageField.H"
#include "LagrangianFields.H"
#include "LagrangianSubFields.H"
#include "volMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class MethodWeightSum, class MethodNoWeightSum, class ... Args>
void Foam::CloudAverageField<Type>::op
(
    MethodWeightSum mws,
    MethodNoWeightSum mnws,
    const LagrangianSubMesh& subMesh,
    const Args& ... args
) const
{
    tcellWeightSum_.valid()
  ? (
        notNull(weightPsiOrPsiState_)
      ? (psiAveragePtr_().*mws)
        (
            subMesh.sub(weightPsiOrPsiState_),
            args ...
        )
     // isNull(weightPsiOrPsiState_)
      : (psiAveragePtr_().*mws)
        (
            toSubField(weightPsiOrPsiDerived_(subMesh)),
            args ...
        )
    )
  : (
        notNull(weightState_) && notNull(weightPsiOrPsiState_)
      ? (psiAveragePtr_().*mnws)
        (
            subMesh.sub(weightState_),
            subMesh.sub(weightPsiOrPsiState_),
            args ...
        )
      : notNull(weightState_) && isNull(weightPsiOrPsiState_)
      ? (psiAveragePtr_().*mnws)
        (
            subMesh.sub(weightState_),
            toSubField(weightPsiOrPsiDerived_(subMesh)),
            args ...
        )
      : isNull(weightState_) && notNull(weightPsiOrPsiState_)
      ? (psiAveragePtr_().*mnws)
        (
            toSubField(weightDerived_(subMesh)),
            subMesh.sub(weightPsiOrPsiState_),
            args ...
        )
     // isNull(weightState_) && isNull(weightPsiOrPsiState_)
      : (psiAveragePtr_().*mnws)
        (
            toSubField(weightDerived_(subMesh)),
            toSubField(weightPsiOrPsiDerived_(subMesh)),
            args ...
        )
    );
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::CloudAverageField<Type>::interpolate
(
    const LagrangianModelRef&,
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianMesh& mesh = subMesh.mesh();

    if (psiAverageState_ == psiAverageState::removed)
    {
        FatalErrorInFunction
            << "Cannot interpolate an average that has removed elements"
            << exit(FatalError);
    }

    if (!psiAveragePtr_.valid())
    {
        psiAveragePtr_ =
            tcellWeightSum_.valid()
          ? (
                notNull(weightPsiOrPsiState_)
              ? LagrangianAverage<Type>::New
                (
                    word(mesh.schemes().averaging(name_)),
                    "average(" + name_ + ')',
                    tcellWeightSum_(),
                    weightPsiOrPsiState_(mesh)
                )
             // isNull(weightPsiOrPsiState_)
              : LagrangianAverage<Type>::New
                (
                    word(mesh.schemes().averaging(name_)),
                    "average(" + name_ + ')',
                    tcellWeightSum_(),
                    weightPsiOrPsiDerived_(mesh)()
                )
            )
          : (
                notNull(weightState_) && notNull(weightPsiOrPsiState_)
              ? LagrangianAverage<Type>::New
                (
                    word(mesh.schemes().averaging(name_)),
                    "average(" + name_ + ')',
                    weightState_(mesh),
                    weightPsiOrPsiState_(mesh)
                )
              : notNull(weightState_) && isNull(weightPsiOrPsiState_)
              ? LagrangianAverage<Type>::New
                (
                    word(mesh.schemes().averaging(name_)),
                    "average(" + name_ + ')',
                    weightState_(mesh),
                    weightPsiOrPsiDerived_(mesh)()
                )
              : isNull(weightState_) && notNull(weightPsiOrPsiState_)
              ? LagrangianAverage<Type>::New
                (
                    word(mesh.schemes().averaging(name_)),
                    "average(" + name_ + ')',
                    weightDerived_(mesh)(),
                    weightPsiOrPsiState_(mesh)
                )
             // isNull(weightState_) && isNull(weightPsiOrPsiState_)
              : LagrangianAverage<Type>::New
                (
                    word(mesh.schemes().averaging(name_)),
                    "average(" + name_ + ')',
                    weightDerived_(mesh)(),
                    weightPsiOrPsiDerived_(mesh)()
                )
            );

        switch (psiAverageState_)
        {
            case psiAverageState::removed:
                op(removeWeightSum_, removeNoWeightSum_, subMesh);
                break;
            case psiAverageState::complete:
                break;
            case psiAverageState::cached:
                op(removeWeightSum_, removeNoWeightSum_, subMesh);
                op(addWeightSum_, addNoWeightSum_, subMesh, true);
                break;
        }
    }

    return psiAveragePtr_->interpolate(subMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CloudAverageField<Type>::CloudAverageField
(
    const word& name,
    const DimensionedField<scalar, volMesh>& cellWeightSum,
    const CloudDerivedField<Type>& weightPsi
)
:
    CloudDerivedField<Type>(*this, &CloudAverageField::interpolate),
    name_(name),
    tcellWeightSum_(cellWeightSum),
    weightState_(NullObjectRef<CloudStateField<scalar>>()),
    weightDerived_(NullObjectRef<CloudDerivedField<scalar>>()),
    weightPsiOrPsiState_(NullObjectRef<CloudStateField<Type>>()),
    weightPsiOrPsiDerived_(weightPsi),
    psiAveragePtr_(nullptr),
    psiAverageState_(psiAverageState::complete)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::CloudAverageField<Type>::remove(const LagrangianSubMesh& subMesh)
{
    psiAverageState_ = psiAverageState::removed;

    if (!psiAveragePtr_.valid()) return;

    CloudDerivedField<Type>::clear(false);

    op(removeWeightSum_, removeNoWeightSum_, subMesh);
}


template<class Type>
void Foam::CloudAverageField<Type>::add
(
    const LagrangianSubMesh& subMesh,
    const bool cache
)
{
    psiAverageState_ =
        cache ? psiAverageState::cached : psiAverageState::complete;

    if (!psiAveragePtr_.valid()) return;

    CloudDerivedField<Type>::clear(false);

    op(addWeightSum_, addNoWeightSum_, subMesh, cache);
}


template<class Type>
void Foam::CloudAverageField<Type>::correct
(
    const LagrangianSubMesh& subMesh,
    const bool cache
)
{
    psiAverageState_ =
        cache ? psiAverageState::cached : psiAverageState::complete;

    if (!psiAveragePtr_.valid()) return;

    CloudDerivedField<Type>::clear(!cache);

    op(correctWeightSum_, correctNoWeightSum_, subMesh, cache);
}


template<class Type>
void Foam::CloudAverageField<Type>::reset()
{
    CloudDerivedField<Type>::clear(true);

    psiAveragePtr_.clear();
    psiAverageState_ = psiAverageState::complete;
}


// ************************************************************************* //

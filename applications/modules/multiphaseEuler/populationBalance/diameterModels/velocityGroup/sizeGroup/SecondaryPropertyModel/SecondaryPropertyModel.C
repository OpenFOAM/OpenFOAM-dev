/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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

#include "SecondaryPropertyModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ModelType>
const Foam::diameterModels::SecondaryPropertyModel<ModelType>&
Foam::diameterModels::SecondaryPropertyModel<ModelType>::correspondingModel
(
    const sizeGroup& fi
) const
{
    const ModelType& model = ModelType::model(fi);

    if (!isA<SecondaryPropertyModel<ModelType>>(model))
    {
        FatalErrorInFunction
            << "Phase " << this->group().phase().name() << " of population "
            << "balance " << this->group().group().popBal().name() << " has a "
            << type() << " " << ModelType::typeName << " but phase "
            << fi.phase().name() << " does not. The " << type() << " "
            << ModelType::typeName << " requires all phases of the population "
            << "balance to " << " have the " << type() << " "
            << ModelType::typeName << "." << exit(FatalError);
    }

    return refCast<const SecondaryPropertyModel<ModelType>>(model);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
Foam::diameterModels::SecondaryPropertyModel<ModelType>::SecondaryPropertyModel
(
    const sizeGroup& group
)
:
    ModelType(group)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::diameterModels::SecondaryPropertyModel<ModelType>::
~SecondaryPropertyModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ModelType>
const Foam::diameterModels::SecondaryPropertyModel<ModelType>&
Foam::diameterModels::SecondaryPropertyModel<ModelType>::model
(
    const sizeGroup& fi
)
{
    const ModelType& model = ModelType::model(fi);

    if (!isA<SecondaryPropertyModel<ModelType>>(model))
    {
        FatalErrorInFunction
            << "Phase " << fi.phase().name() << " does not have a "
            << ModelType::typeName << " with a secondary property"
            << exit(FatalError);
    }

    return refCast<const SecondaryPropertyModel<ModelType>>(model);
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::addCoalescence
(
    const volScalarField::Internal& Su,
    const sizeGroup& fj,
    const sizeGroup& fk
)
{
    const volScalarField::Internal& propj = correspondingModel(fj).fld();
    const volScalarField::Internal& propk = correspondingModel(fk).fld();

    src() += (propj*fj.x() + propk*fk.x())/(fj.x() + fk.x())*Su;
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::addBreakup
(
    const volScalarField::Internal& Su,
    const sizeGroup& fj
)
{
    const volScalarField::Internal& propj = correspondingModel(fj).fld();

    src() += propj*Su;
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::reset()
{
    src() = Zero;
}


// ************************************************************************* //

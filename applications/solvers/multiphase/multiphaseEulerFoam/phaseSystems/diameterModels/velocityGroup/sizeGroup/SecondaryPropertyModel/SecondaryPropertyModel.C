/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
#include "shapeModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
Foam::diameterModels::SecondaryPropertyModel<ModelType>::SecondaryPropertyModel
(
    const dictionary& dict,
    const sizeGroup& group
)
:
    ModelType(dict, group),
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeid(ModelType).name(), group.name()),
            group.time().constant(),
            group.mesh()
        )
    ),
    sizeGroup_(group),
    SecondaryPropertyModelTable_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::diameterModels::SecondaryPropertyModel<ModelType>::
~SecondaryPropertyModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ModelType>
const typename Foam::diameterModels::SecondaryPropertyModel<ModelType>::SpTable&
Foam::diameterModels::SecondaryPropertyModel<ModelType>::
SecondaryPropertyModelTable()
{
    if (SecondaryPropertyModelTable_.empty())
    {
        SecondaryPropertyModelTable_ =
            sizeGroup_.mesh().template lookupClass
            <
                SecondaryPropertyModel<ModelType>
            >();
    }

    return SecondaryPropertyModelTable_;
}


template<class ModelType>
const Foam::word Foam::diameterModels::SecondaryPropertyModel<ModelType>::
SecondaryPropertyName
(
    const sizeGroup& fi
) const
{
    return word(IOobject::groupName(typeid(ModelType).name(), fi.name()));
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::addCoalescence
(
    const volScalarField& Su,
    const sizeGroup& fj,
    const sizeGroup& fk
)
{
    const volScalarField& propj =
        SecondaryPropertyModelTable()[SecondaryPropertyName(fj)]->fld();

    const volScalarField& propk =
        SecondaryPropertyModelTable()[SecondaryPropertyName(fk)]->fld();

    src() += (propj*fj.x() + propk*fk.x())/(fj.x() + fk.x())*Su;
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::addBreakup
(
    const volScalarField& Su,
    const sizeGroup& fj
)
{
    const volScalarField& propj =
        SecondaryPropertyModelTable()[SecondaryPropertyName(fj)]->fld();

    src() += propj*Su;
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::addDrift
(
    const volScalarField& Su,
    const sizeGroup& fu,
    const driftModel& model
)
{
    const volScalarField& propu =
        SecondaryPropertyModelTable()[SecondaryPropertyName(fu)]->fld();

    src() += propu*Su;
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::addNucleation
(
    const volScalarField& Su,
    const sizeGroup& fi,
    const nucleationModel& model
)
{
    const volScalarField& prop =
        SecondaryPropertyModelTable()[SecondaryPropertyName(fi)]->fld();

    src() += prop*Su;
}


template<class ModelType>
void Foam::diameterModels::SecondaryPropertyModel<ModelType>::reset()
{
    src() = Zero;
}


template<class ModelType>
bool Foam::diameterModels::SecondaryPropertyModel<ModelType>::
writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //

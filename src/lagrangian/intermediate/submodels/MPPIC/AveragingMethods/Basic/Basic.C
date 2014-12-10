/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "Basic.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Basic<Type>::Basic
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    AveragingMethod<Type>(io, dict, mesh, labelList(1, mesh.nCells())),
    data_(FieldField<Field, Type>::operator[](0)),
    dataGrad_(mesh.nCells())
{}


template<class Type>
Foam::AveragingMethods::Basic<Type>::Basic
(
    const Basic<Type>& am
)
:
    AveragingMethod<Type>(am),
    data_(FieldField<Field, Type>::operator[](0)),
    dataGrad_(am.dataGrad_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Basic<Type>::~Basic()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Basic<Type>::updateGrad()
{
    GeometricField<Type, fvPatchField, volMesh> tempData
    (
        IOobject
        (
            "BasicAverage::Data",
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero),
        zeroGradientFvPatchField<Type>::typeName
    );
    tempData.internalField() = data_;
    tempData.correctBoundaryConditions();
    dataGrad_ = fvc::grad(tempData)->internalField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Basic<Type>::add
(
    const point position,
    const tetIndices& tetIs,
    const Type& value
)
{
    data_[tetIs.cell()] += value/this->mesh_.V()[tetIs.cell()];
}


template<class Type>
Type Foam::AveragingMethods::Basic<Type>::interpolate
(
    const point position,
    const tetIndices& tetIs
) const
{
    return data_[tetIs.cell()];
}


template<class Type>
typename Foam::AveragingMethods::Basic<Type>::TypeGrad
Foam::AveragingMethods::Basic<Type>::interpolateGrad
(
    const point position,
    const tetIndices& tetIs
) const
{
    return dataGrad_[tetIs.cell()];
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AveragingMethods::Basic<Type>::internalField() const
{
    return tmp<Field<Type> >(data_);
}


// ************************************************************************* //

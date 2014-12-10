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

#include "Moment.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Moment<Type>::Moment
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    AveragingMethod<Type>(io, dict, mesh, labelList(4, mesh.nCells())),
    data_(FieldField<Field, Type>::operator[](0)),
    dataX_(FieldField<Field, Type>::operator[](1)),
    dataY_(FieldField<Field, Type>::operator[](2)),
    dataZ_(FieldField<Field, Type>::operator[](3)),
    transform_(mesh.nCells(), symmTensor::zero),
    scale_(0.5*pow(mesh.V(), 1.0/3.0))
{
    scalar a = 1.0/24.0;
    scalar b = 0.5854101966249685;
    scalar c = 0.1381966011250105;

    scalarField wQ(4);
    wQ[0] = a;
    wQ[1] = a;
    wQ[2] = a;
    wQ[3] = a;

    vectorField xQ(4);
    xQ[0] = vector(b, c, c);
    xQ[1] = vector(c, b, c);
    xQ[2] = vector(c, c, b);
    xQ[3] = vector(c, c, c);

    forAll(mesh.C(), cellI)
    {
        const List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(mesh, cellI);

        symmTensor A(symmTensor::zero);

        forAll(cellTets, tetI)
        {
            const tetIndices& tetIs = cellTets[tetI];
            const label faceI = tetIs.face();
            const face& f = mesh.faces()[faceI];

            const tensor T
            (
                tensor
                (
                    mesh.points()[f[tetIs.faceBasePt()]] - mesh.C()[cellI],
                    mesh.points()[f[tetIs.facePtA()]] - mesh.C()[cellI],
                    mesh.points()[f[tetIs.facePtB()]] - mesh.C()[cellI]
                ).T()
            );

            const vectorField d((T & xQ)/scale_[cellI]);

            const scalar v(6.0*tetIs.tet(mesh).mag()/mesh.V()[cellI]);

            A += v*sum(wQ*sqr(d));
        }

        transform_[cellI] = inv(A);
    }
}


template<class Type>
Foam::AveragingMethods::Moment<Type>::Moment
(
    const Moment<Type>& am
)
:
    AveragingMethod<Type>(am),
    data_(FieldField<Field, Type>::operator[](0)),
    dataX_(FieldField<Field, Type>::operator[](1)),
    dataY_(FieldField<Field, Type>::operator[](2)),
    dataZ_(FieldField<Field, Type>::operator[](3)),
    transform_(am.transform_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Moment<Type>::~Moment()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Moment<Type>::updateGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Moment<Type>::add
(
    const point position,
    const tetIndices& tetIs,
    const Type& value
)
{
    const label cellI = tetIs.cell();

    const Type v = value/this->mesh_.V()[cellI];
    const TypeGrad dv =
        transform_[cellI]
      & (
            v
          * (position - this->mesh_.C()[cellI])
          / scale_[cellI]
        );

    data_[cellI] += v;
    dataX_[cellI] += v + dv.x();
    dataY_[cellI] += v + dv.y();
    dataZ_[cellI] += v + dv.z();
}


template<class Type>
Type Foam::AveragingMethods::Moment<Type>::interpolate
(
    const point position,
    const tetIndices& tetIs
) const
{
    const label cellI = tetIs.cell();

    return
        data_[cellI]
      + (
            TypeGrad
            (
                dataX_[cellI] - data_[cellI],
                dataY_[cellI] - data_[cellI],
                dataZ_[cellI] - data_[cellI]
            )
          & (position - this->mesh_.C()[cellI])
          / scale_[cellI]
        );
}


template<class Type>
typename Foam::AveragingMethods::Moment<Type>::TypeGrad
Foam::AveragingMethods::Moment<Type>::interpolateGrad
(
    const point position,
    const tetIndices& tetIs
) const
{
    const label cellI(tetIs.cell());

    return
        TypeGrad
        (
            dataX_[cellI] - data_[cellI],
            dataY_[cellI] - data_[cellI],
            dataZ_[cellI] - data_[cellI]
        )/scale_[cellI];
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AveragingMethods::Moment<Type>::internalField() const
{
    return tmp<Field<Type> >(data_);
}


// ************************************************************************* //

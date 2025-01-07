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

#include "cellPoint_LagrangianAverage.H"
#include "cellPointLagrangianAddressor.H"
#include "LagrangianFields.H"
#include "oneField.H"
#include "tetIndices.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class CellWeight>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::LagrangianAverages::cellPoint<Type>::cellVolumeWeightedSum
(
    const polyMesh& pMesh,
    const CellWeight& cellWeight
)
{
    tmp<Field<scalar>> tresult(new Field<scalar>(pMesh.nCells(), scalar(0)));
    Field<scalar>& result = tresult.ref();

    forAll(pMesh.cells(), celli)
    {
        const Foam::cell& c = pMesh.cells()[celli];

        forAll(c, cFacei)
        {
            const label facei = c[cFacei];
            const Foam::face& f = pMesh.faces()[facei];

            for (label fTrii = 0; fTrii < f.nTriangles(); ++ fTrii)
            {
                const tetIndices tetPoints(celli, facei, fTrii + 1);
                const scalar v = tetPoints.tet(pMesh).mag();

                result[celli] += v*cellWeight[celli];
            }
        }
    }

    return tresult;
}


template<class Type>
template<class CellWeight>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::LagrangianAverages::cellPoint<Type>::pointVolumeWeightedSum
(
    const polyMesh& pMesh,
    const CellWeight& cellWeight
)
{
    tmp<Field<scalar>> tresult(new Field<scalar>(pMesh.nPoints(), scalar(0)));
    Field<scalar>& result = tresult.ref();

    forAll(pMesh.cells(), celli)
    {
        const Foam::cell& c = pMesh.cells()[celli];

        forAll(c, cFacei)
        {
            const label facei = c[cFacei];
            const Foam::face& f = pMesh.faces()[facei];

            for (label fTrii = 0; fTrii < f.nTriangles(); ++ fTrii)
            {
                const tetIndices tetPoints(celli, facei, fTrii + 1);
                const scalar v = tetPoints.tet(pMesh).mag();
                const triFace triPoints = tetPoints.faceTriIs(pMesh);

                forAll(triPoints, triPointi)
                {
                    const label pointi = triPoints[triPointi];

                    result[pointi] += v*cellWeight[celli];
                }
            }
        }
    }

    syncTools::syncPointList
    (
        pMesh,
        result,
        plusEqOp<scalar>(),
        scalar(0)
    );

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::LagrangianAverages::cellPoint<Type>::cellWeightSum
(
    const polyMesh& pMesh,
    const Field<scalar>& weightSum
)
{
    return weightSum/4;
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::LagrangianAverages::cellPoint<Type>::pointWeightSum
(
    const polyMesh& pMesh,
    const Field<scalar>& weightSum
)
{
    return
        pointVolumeWeightedSum
        (
            pMesh,
            (weightSum/cellVolumeWeightedSum(pMesh, oneField()))()
        )/4;
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::clear(data& d)
{
    forAll(d.cellAvgCell_, cellAvgi)
    {
        d.cellCellAvg_[d.cellAvgCell_[cellAvgi]] = -1;
    }
    forAll(d.pointAvgPoint_, pointAvgi)
    {
        d.pointPointAvg_[d.pointAvgPoint_[pointAvgi]] = -1;
    }
    d.cellAvgCell_.clear();
    d.pointAvgPoint_.clear();
    d.cellAvgCount_.clear();
    d.pointAvgCount_.clear();
    if (d.cellAvgWeightSumPtr_.valid()) d.cellAvgWeightSumPtr_->clear();
    if (d.pointAvgWeightSumPtr_.valid()) d.pointAvgWeightSumPtr_->clear();
    d.cellAvgSum_.clear();
    d.pointAvgSum_.clear();
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::removeFromCells
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi,
    data& d
)
{
    const LagrangianSubMesh& subMesh = psiOrWeightPsi.mesh();
    const LagrangianMesh& mesh = subMesh.mesh();

    if (notNull(weightOrNull) != d.cellAvgWeightSumPtr_.valid())
    {
        FatalErrorInFunction
            << "Inconsistent weight specification for average of field "
            << psiOrWeightPsi.name() << exit(FatalError);
    }

    // Remove from the cell averages
    forAll(psiOrWeightPsi, subi)
    {
        const label i = subMesh.start() + subi;

        const barycentric& coordinates = mesh.coordinates()[i];
        const label celli = mesh.celli()[i];

        const label cellAvgi = d.cellCellAvg_[celli];

        // Check that this cell can be removed from
        if (cellAvgi == -1 || d.cellAvgCount_[cellAvgi] == 0)
        {
            FatalErrorInFunction
                << "Negative cell count for average of field "
                << psiOrWeightPsi.name() << exit(FatalError);
        }

        // Remove from the cell average
        d.cellAvgCount_[cellAvgi] --;
        if (notNull(weightOrNull))
        {
            d.cellAvgWeightSumPtr_()[cellAvgi] -=
                coordinates.a()*weightOrNull[subi];
            d.cellAvgSum_[cellAvgi] -=
                coordinates.a()*weightOrNull[subi]*psiOrWeightPsi[subi];
        }
        else
        {
            d.cellAvgSum_[cellAvgi] -= coordinates.a()*psiOrWeightPsi[subi];
        }
    }
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::addToCells
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi,
    data& d
)
{
    const LagrangianSubMesh& subMesh = psiOrWeightPsi.mesh();
    const LagrangianMesh& mesh = subMesh.mesh();

    if (notNull(weightOrNull) != d.cellAvgWeightSumPtr_.valid())
    {
        FatalErrorInFunction
            << "Inconsistent weight specification for average of field "
            << psiOrWeightPsi.name() << exit(FatalError);
    }

    // Add to the cell averages
    forAll(psiOrWeightPsi, subi)
    {
        const label i = subMesh.start() + subi;

        const barycentric& coordinates = mesh.coordinates()[i];
        const label celli = mesh.celli()[i];

        label cellAvgi = d.cellCellAvg_[celli];

        // Initialise this cell if it is newly a part of the average
        if (cellAvgi == -1)
        {
            cellAvgi = d.cellAvgCell_.size();
            d.cellCellAvg_[celli] = cellAvgi;
            d.cellAvgCell_.append(celli);
            d.cellAvgCount_.append(label(0));
            if (notNull(weightOrNull))
            {
                d.cellAvgWeightSumPtr_().append(scalar(0));
            }
            d.cellAvgSum_.append(pTraits<Type>::zero);
        }

        // Add to the cell average
        d.cellAvgCount_[cellAvgi] ++;
        if (notNull(weightOrNull))
        {
            d.cellAvgWeightSumPtr_()[cellAvgi] +=
                coordinates.a()*weightOrNull[subi];
            d.cellAvgSum_[cellAvgi] +=
                coordinates.a()*weightOrNull[subi]*psiOrWeightPsi[subi];
        }
        else
        {
            d.cellAvgSum_[cellAvgi] +=
                coordinates.a()*psiOrWeightPsi[subi];
        }
    }
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::addToPoints
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi,
    data& d
)
{
    const LagrangianSubMesh& subMesh = psiOrWeightPsi.mesh();
    const LagrangianMesh& mesh = subMesh.mesh();

    if (notNull(weightOrNull) != d.cellAvgWeightSumPtr_.valid())
    {
        FatalErrorInFunction
            << "Inconsistent weight specification for average of field "
            << psiOrWeightPsi.name() << exit(FatalError);
    }

    // Add to the point averages
    forAll(psiOrWeightPsi, subi)
    {
        const label i = subMesh.start() + subi;

        const barycentric& coordinates = mesh.coordinates()[i];
        const label celli = mesh.celli()[i];
        const label facei = mesh.facei()[i];
        const label faceTrii = mesh.faceTrii()[i];

        const triFace triPoints =
            tetIndices(celli, facei, faceTrii).faceTriIs(mesh.mesh());

        forAll(triPoints, triPointi)
        {
            const label pointi = triPoints[triPointi];

            label pointAvgi = d.pointPointAvg_[pointi];

            // Initialise this point if it is newly a part of the average
            if (pointAvgi == -1)
            {
                pointAvgi = d.pointAvgPoint_.size();
                d.pointPointAvg_[pointi] = pointAvgi;
                d.pointAvgPoint_.append(pointi);
                d.pointAvgCount_.append(label(0));
                if (notNull(weightOrNull))
                {
                    d.pointAvgWeightSumPtr_().append(scalar(0));
                }
                d.pointAvgSum_.append(pTraits<Type>::zero);
            }

            // Add to the point average
            d.pointAvgCount_[pointAvgi] ++;
            if (notNull(weightOrNull))
            {
                d.pointAvgWeightSumPtr_()[pointAvgi] +=
                    coordinates[triPointi + 1]*weightOrNull[subi];
                d.pointAvgSum_[pointAvgi] +=
                    coordinates[triPointi + 1]
                   *weightOrNull[subi]
                   *psiOrWeightPsi[subi];
            }
            else
            {
                d.pointAvgSum_[pointAvgi] +=
                    coordinates[triPointi + 1]*psiOrWeightPsi[subi];
            }
        }
    }

    // Expand to include any points that are being added to remotely
    cellPointLagrangianAddressor::New(mesh.mesh()).sync
    (
        d.pointPointAvg_,
        d.pointAvgPoint_
    );
    d.pointAvgCount_.resize(d.pointAvgPoint_.size(), label(0));
    if (notNull(weightOrNull))
    {
        d.pointAvgWeightSumPtr_().resize(d.pointAvgPoint_.size(), scalar(0));
    }
    d.pointAvgSum_.resize(d.pointAvgPoint_.size(), pTraits<Type>::zero);

    // Sum on coupled points
    syncTools::syncPointList
    (
        mesh.mesh(),
        d.pointAvgPoint_,
        d.pointAvgCount_,
        plusEqOp<label>(),
        label(0)
    );
    if (notNull(weightOrNull))
    {
        syncTools::syncPointList
        (
            mesh.mesh(),
            d.pointAvgPoint_,
            d.pointAvgWeightSumPtr_(),
            plusEqOp<scalar>(),
            scalar(0)
        );
    }
    syncTools::syncPointList
    (
        mesh.mesh(),
        d.pointAvgPoint_,
        d.pointAvgSum_,
        plusEqOp<Type>(),
        pTraits<Type>::zero
    );
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::removeFromPoints
(
    const data& dd,
    data& d
)
{
    // Remove from the point averages
    forAll(dd.pointAvgPoint_, dPointAvgi)
    {
        const label pointi = dd.pointAvgPoint_[dPointAvgi];
        const label pointAvgi = d.pointPointAvg_[pointi];

        // Check that this point can be removed from
        if
        (
            pointAvgi == -1
         || d.pointAvgCount_[pointAvgi] < dd.pointAvgCount_[dPointAvgi]
        )
        {
            FatalErrorInFunction
                << "Negative point count for average "
                << exit(FatalError);
        }

        // Remove from the point average
        d.pointAvgCount_[pointAvgi] -= dd.pointAvgCount_[dPointAvgi];
        if (d.pointAvgWeightSumPtr_.valid())
        {
            d.pointAvgWeightSumPtr_()[pointAvgi] -=
                dd.pointAvgWeightSumPtr_()[dPointAvgi];
        }
        d.pointAvgSum_[pointAvgi] -= dd.pointAvgSum_[dPointAvgi];
    }
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::addToPoints
(
    const data& dd,
    data& d
)
{
    forAll(dd.pointAvgPoint_, dPointAvgi)
    {
        const label pointi = dd.pointAvgPoint_[dPointAvgi];

        label pointAvgi = d.pointPointAvg_[pointi];

        // Initialise this point if it is newly a part of the average
        if (pointAvgi == -1)
        {
            pointAvgi = d.pointAvgPoint_.size();
            d.pointPointAvg_[pointi] = pointAvgi;
            d.pointAvgPoint_.append(pointi);
            d.pointAvgCount_.append(label(0));
            if (d.pointAvgWeightSumPtr_.valid())
            {
                d.pointAvgWeightSumPtr_().append(scalar(0));
            }
            d.pointAvgSum_.append(pTraits<Type>::zero);
        }

        // Add to the point average
        d.pointAvgCount_[pointAvgi] += dd.pointAvgCount_[dPointAvgi];
        if (d.pointAvgWeightSumPtr_.valid())
        {
            d.pointAvgWeightSumPtr_()[pointAvgi] +=
                dd.pointAvgWeightSumPtr_()[dPointAvgi];
        }
        d.pointAvgSum_[pointAvgi] += dd.pointAvgSum_[dPointAvgi];
    }
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::interpolate
(
    LagrangianSubField<Type>& result
) const
{
    const LagrangianSubMesh& subMesh = result.mesh();
    const LagrangianMesh& mesh = subMesh.mesh();

    forAll(subMesh, subi)
    {
        const label i = subMesh.start() + subi;

        const barycentric& coordinates = mesh.coordinates()[i];
        const label celli = mesh.celli()[i];
        const label facei = mesh.facei()[i];
        const label faceTrii = mesh.faceTrii()[i];

        const triFace triPoints =
            tetIndices(celli, facei, faceTrii).faceTriIs(mesh.mesh());

        const label cellAvgi = data_.cellCellAvg_[celli];
        const label dCellAvgi = dData_.cellCellAvg_[celli];
        const bool haveData =
            cellAvgi != -1 && data_.cellAvgCount_[cellAvgi] != 0;
        const bool haveDData =
            dCellAvgi != -1 && dData_.cellAvgCount_[dCellAvgi] != 0;

        if (!haveData && !haveDData)
        {
            FatalErrorInFunction
                << "Interpolated value requested for a cell in which no "
                << "elements have been averaged"
                << exit(FatalError);
        }

        static const Type& z = pTraits<Type>::zero;

        Type wr =
            coordinates.a()
           *(
                (haveData ? data_.cellAvgSum_[cellAvgi] : z)
              + (haveDData ? dData_.cellAvgSum_[dCellAvgi] : z)
            );

       scalar w =
            coordinates.a()
           *(
                !cellWeightSumPtr_.valid()
              ? (haveData ? data_.cellAvgWeightSumPtr_()[cellAvgi] : 0)
              + (haveDData ? dData_.cellAvgWeightSumPtr_()[dCellAvgi] : 0)
              : cellWeightSumPtr_()[celli]
            );

        forAll(triPoints, triPointi)
        {
            const label pointi = triPoints[triPointi];

            const label pointAvgi = data_.pointPointAvg_[pointi];
            const label dPointAvgi = dData_.pointPointAvg_[pointi];
            const bool haveData =
                pointAvgi != -1 && data_.pointAvgCount_[pointAvgi] != 0;
            const bool haveDData =
                dPointAvgi != -1 && dData_.pointAvgCount_[dPointAvgi] != 0;

            wr +=
                coordinates[triPointi + 1]
               *(
                    (haveData ? data_.pointAvgSum_[pointAvgi] : z)
                  + (haveDData ? dData_.pointAvgSum_[dPointAvgi] : z)
                );

            w +=
                coordinates[triPointi + 1]
               *(
                    !pointWeightSumPtr_.valid()
                  ? (haveData ? data_.pointAvgWeightSumPtr_()[pointAvgi] : 0)
                  + (haveDData ? dData_.pointAvgWeightSumPtr_()[dPointAvgi] : 0)
                  : pointWeightSumPtr_()[pointi]
                );
        }

        result[subi] = wr/w;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianAverages::cellPoint<Type>::data::data
(
    const label nCells,
    const label nPoints,
    const bool hasWeightSum
)
:
    cellCellAvg_(nCells, label(-1)),
    pointPointAvg_(nPoints, label(-1)),
    cellAvgCell_(),
    pointAvgPoint_(),
    cellAvgCount_(),
    pointAvgCount_(),
    cellAvgWeightSumPtr_(hasWeightSum ? new DynamicList<scalar>() : nullptr),
    pointAvgWeightSumPtr_(hasWeightSum ? new DynamicList<scalar>() : nullptr),
    cellAvgSum_(),
    pointAvgSum_()
{}


template<class Type>
Foam::LagrangianAverages::cellPoint<Type>::cellPoint
(
    const word& name,
    const LagrangianMesh& mesh,
    const dimensionSet& dimensions,
    const Field<scalar>& weightSum
)
:
    LagrangianAverage<Type>(name, mesh, dimensions),
    cellWeightSumPtr_
    (
        !isNull(weightSum)
      ? cellWeightSum(mesh.mesh(), weightSum).ptr()
      : nullptr
    ),
    pointWeightSumPtr_
    (
        !isNull(weightSum)
      ? pointWeightSum(mesh.mesh(), weightSum).ptr()
      : nullptr
    ),
    data_(mesh.mesh().nCells(), mesh.mesh().nPoints(), isNull(weightSum)),
    dDataIsValid_(false),
    dData_(mesh.mesh().nCells(), mesh.mesh().nPoints(), isNull(weightSum))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianAverages::cellPoint<Type>::~cellPoint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::remove
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi
)
{
    if (dDataIsValid_)
    {
        FatalErrorInFunction
            << "Cannot remove from an average with a cached difference"
            << exit(FatalError);
    }

    removeFromCells(weightOrNull, psiOrWeightPsi, data_);

    addToPoints(weightOrNull, psiOrWeightPsi, dData_);

    removeFromPoints(dData_, data_);

    clear(dData_);
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::add
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi,
    const bool cache
)
{
    if (dDataIsValid_)
    {
        FatalErrorInFunction
            << "Cannot add to an average with a cached difference"
            << exit(FatalError);
    }

    dDataIsValid_ = cache;

    addToCells(weightOrNull, psiOrWeightPsi, cache ? dData_ : data_);

    addToPoints(weightOrNull, psiOrWeightPsi, dData_);

    if (!cache)
    {
        addToPoints(dData_, data_);

        clear(dData_);
    }
}


template<class Type>
void Foam::LagrangianAverages::cellPoint<Type>::correct
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi,
    const bool cache
)
{
    if (!dDataIsValid_)
    {
        FatalErrorInFunction
            << "Cannot correct an average without a cached difference"
            << exit(FatalError);
    }

    clear(dData_);

    dDataIsValid_ = cache;

    addToCells(weightOrNull, psiOrWeightPsi, cache ? dData_ : data_);

    addToPoints(weightOrNull, psiOrWeightPsi, dData_);

    if (!cache)
    {
        addToPoints(dData_, data_);

        clear(dData_);
    }
}


// ************************************************************************* //

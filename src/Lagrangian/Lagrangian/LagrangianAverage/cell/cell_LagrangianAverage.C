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

#include "cell_LagrangianAverage.H"
#include "LagrangianFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::LagrangianAverages::cell<Type>::clear(data& d)
{
    forAll(d.cellAvgCell_, cellAvgi)
    {
        d.cellCellAvg_[d.cellAvgCell_[cellAvgi]] = -1;
    }
    d.cellAvgCell_.clear();
    d.cellAvgCount_.clear();
    if (d.cellAvgWeightSumPtr_.valid()) d.cellAvgWeightSumPtr_->clear();
    d.cellAvgSum_.clear();
}


template<class Type>
void Foam::LagrangianAverages::cell<Type>::remove
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi,
    data& d
)
{
    const LagrangianSubMesh& subMesh = psiOrWeightPsi.mesh();

    //Info<< "*** Removing sub mesh with " << subMesh.size()
    //    << " elements" << endl;

    if (notNull(weightOrNull) != d.cellAvgWeightSumPtr_.valid())
    {
        FatalErrorInFunction
            << "Inconsistent weight specification for average of field "
            << psiOrWeightPsi.name() << exit(FatalError);
    }

    forAll(psiOrWeightPsi, subi)
    {
        const label celli = subMesh.mesh().celli()[subMesh.start() + subi];
        const label cellAvgi = d.cellCellAvg_[celli];

        //Info<< " -> Remove from cell #" << celli << ", leaving "
        //    << d.cellAvgCount_[cellAvgi] - 1 << " samples" << endl;

        if (cellAvgi == -1 || d.cellAvgCount_[cellAvgi] == 0)
        {
            FatalErrorInFunction
                << "Negative count for average of field "
                << psiOrWeightPsi.name() << exit(FatalError);
        }

        // Remove from the cell average
        d.cellAvgCount_[cellAvgi] --;
        if (notNull(weightOrNull))
        {
            d.cellAvgWeightSumPtr_()[cellAvgi] -= weightOrNull[subi];
            d.cellAvgSum_[cellAvgi] -= weightOrNull[subi]*psiOrWeightPsi[subi];
        }
        else
        {
            d.cellAvgSum_[cellAvgi] -= psiOrWeightPsi[subi];
        }
    }

    //forAll(dcellAvgCell__, cellAvgi)
    //{
    //    Info<< " !! Cell #" << d.cellAvgCell_[cellAvgi] << " has "
    //        << d.cellAvgCount_[cellAvgi] << " samples" << endl;
    //}
}


template<class Type>
void Foam::LagrangianAverages::cell<Type>::add
(
    const LagrangianSubSubField<scalar>& weightOrNull,
    const LagrangianSubSubField<Type>& psiOrWeightPsi,
    data& d
)
{
    const LagrangianSubMesh& subMesh = psiOrWeightPsi.mesh();

    //Info<< "*** Adding sub mesh with " << subMesh.size()
    //    << " elements" << endl;

    if (notNull(weightOrNull) != d.cellAvgWeightSumPtr_.valid())
    {
        FatalErrorInFunction
            << "Inconsistent weight specification for average of field "
            << psiOrWeightPsi.name() << exit(FatalError);
    }

    forAll(psiOrWeightPsi, subi)
    {
        const label celli = subMesh.mesh().celli()[subMesh.start() + subi];
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

        //Info<< " -> Add to cell #" << celli << ", giving "
        //    << d.cellAvgCount_[cellAvgi] + 1 << " samples" << endl;

        // Add to the cell average
        d.cellAvgCount_[cellAvgi] ++;
        if (notNull(weightOrNull))
        {
            d.cellAvgWeightSumPtr_()[cellAvgi] += weightOrNull[subi];
            d.cellAvgSum_[cellAvgi] += weightOrNull[subi]*psiOrWeightPsi[subi];
        }
        else
        {
            d.cellAvgSum_[cellAvgi] += psiOrWeightPsi[subi];
        }
    }

    //forAll(d.cellAvgCell_, cellAvgi)
    //{
    //    Info<< " !! Cell #" << d.cellAvgCell_[cellAvgi] << " has "
    //        << d.cellAvgCount_[cellAvgi] << " samples" << endl;
    //}
}


template<class Type>
void Foam::LagrangianAverages::cell<Type>::interpolate
(
    LagrangianSubField<Type>& result
) const
{
    const LagrangianSubMesh& subMesh = result.mesh();

    forAll(subMesh, subi)
    {
        const label celli = subMesh.mesh().celli()[subMesh.start() + subi];
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

        result[subi] =
            (
                (haveData ? data_.cellAvgSum_[cellAvgi] : z)
              + (haveDData ? dData_.cellAvgSum_[dCellAvgi] : z)
            )
           /(
                isNull(weightSum_)
              ? (haveData ? data_.cellAvgWeightSumPtr_()[cellAvgi] : 0)
              + (haveDData ? dData_.cellAvgWeightSumPtr_()[dCellAvgi] : 0)
              : weightSum_[celli]
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianAverages::cell<Type>::data::data
(
    const label nCells,
    const bool hasWeightSum
)
:
    cellCellAvg_(nCells, label(-1)),
    cellAvgCell_(),
    cellAvgCount_(),
    cellAvgWeightSumPtr_(hasWeightSum ? new DynamicList<scalar>() : nullptr),
    cellAvgSum_()
{}


template<class Type>
Foam::LagrangianAverages::cell<Type>::cell
(
    const word& name,
    const LagrangianMesh& mesh,
    const dimensionSet& dimensions,
    const Field<scalar>& weightSum
)
:
    LagrangianAverage<Type>(name, mesh, dimensions),
    weightSum_(weightSum),
    data_(mesh.mesh().nCells(), isNull(weightSum)),
    dDataIsValid_(false),
    dData_(mesh.mesh().nCells(), isNull(weightSum))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianAverages::cell<Type>::~cell()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::LagrangianAverages::cell<Type>::remove
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

    remove(weightOrNull, psiOrWeightPsi, data_);
}


template<class Type>
void Foam::LagrangianAverages::cell<Type>::add
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

    add(weightOrNull, psiOrWeightPsi, cache ? dData_ : data_);
}


template<class Type>
void Foam::LagrangianAverages::cell<Type>::correct
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

    add(weightOrNull, psiOrWeightPsi, cache ? dData_ : data_);
}


// ************************************************************************* //

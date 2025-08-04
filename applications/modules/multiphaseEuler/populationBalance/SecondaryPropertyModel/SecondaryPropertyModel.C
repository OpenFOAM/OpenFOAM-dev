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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
Foam::populationBalance::SecondaryPropertyModel<ModelType>::
SecondaryPropertyModel
(
    const populationBalanceModel& popBal
)
:
    ModelType(popBal)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::populationBalance::SecondaryPropertyModel<ModelType>::
~SecondaryPropertyModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ModelType>
void Foam::populationBalance::SecondaryPropertyModel<ModelType>::addCoalescence
(
    const volScalarField::Internal& Su,
    const label i,
    const label j,
    const label k
)
{
    const PtrList<dimensionedScalar>& vs = this->popBal().vs();

    const volScalarField::Internal& propj = fld(j);
    const volScalarField::Internal& propk = fld(k);

    src(i) += (propj*vs[j] + propk*vs[k])/(vs[j] + vs[k])*Su;
}


template<class ModelType>
void Foam::populationBalance::SecondaryPropertyModel<ModelType>::addBreakup
(
    const volScalarField::Internal& Su,
    const label i,
    const label j
)
{
    const volScalarField::Internal& propj = fld(j);

    src(i) += propj*Su;
}


template<class ModelType>
void Foam::populationBalance::SecondaryPropertyModel<ModelType>::reset()
{
    forAll(this->popBal().fs(), i)
    {
        src(i) = Zero;
    }
}


// ************************************************************************* //

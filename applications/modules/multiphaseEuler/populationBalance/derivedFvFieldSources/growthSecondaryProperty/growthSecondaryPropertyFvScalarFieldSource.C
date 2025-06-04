/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "growthSecondaryPropertyFvScalarFieldSource.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::growthSecondaryPropertyFvScalarFieldSource::
growthSecondaryPropertyFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    growthFvScalarFieldSource(iF, dict),
    groupPropertyFvScalarField(iF)
{}


Foam::growthSecondaryPropertyFvScalarFieldSource::
growthSecondaryPropertyFvScalarFieldSource
(
    const growthSecondaryPropertyFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    growthFvScalarFieldSource(field, iF),
    groupPropertyFvScalarField(iF)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthSecondaryPropertyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();
    const volScalarField& fi = popBal.f(i);

    return fi.sources()[model.name()].internalCoeff(model, source)*fi;
}


Foam::Pair<Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>>
Foam::growthSecondaryPropertyFvScalarFieldSource::sourceCoeffs
(
    const fvSource& model
) const
{
    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();
    const volScalarField& fi = popBal.f(i);

    Pair<tmp<DimensionedField<scalar, volMesh>>> fiSourceCoeffs =
        refCast<const growthFvScalarFieldSource>
        (
            fi.sources()[model.name()]
        ).sourceCoeffs(model);

    Pair<tmp<DimensionedField<scalar, volMesh>>> tsourceCoeffs;

    if (i != 0)
    {
        tsourceCoeffs.first() = fiSourceCoeffs.first()*value(i - 1, model);
    }

    if (i != popBal.nGroups() - 1)
    {
        tsourceCoeffs.second() = fiSourceCoeffs.second()*value(i + 1, model);
    }

    return tsourceCoeffs;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthSecondaryPropertyFvScalarFieldSource::sourceCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();

    Pair<tmp<DimensionedField<scalar, volMesh>>> tsourceCoeffs =
        sourceCoeffs(model);

    return
        i ==  popBal.diameters()[i].iFirst()
      ? neg(source)*tsourceCoeffs.second()
      : i == popBal.diameters()[i].iLast()
      ? pos(source)*tsourceCoeffs.first()
      : pos(source)*tsourceCoeffs.first()
      + neg(source)*tsourceCoeffs.second();
}


// ************************************************************************* //

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
    secondaryPropertyFvScalarFieldSource(iF)
{}


Foam::growthSecondaryPropertyFvScalarFieldSource::
growthSecondaryPropertyFvScalarFieldSource
(
    const growthSecondaryPropertyFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    growthFvScalarFieldSource(field, iF),
    secondaryPropertyFvScalarFieldSource(iF)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthSecondaryPropertyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const diameterModels::sizeGroup& fi = this->fi();

    return fi.sources()[model.name()].internalCoeff(model, source)*fi;
}


Foam::Pair<Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>>
Foam::growthSecondaryPropertyFvScalarFieldSource::sourceCoeffs
(
    const fvSource& model
) const
{
    const diameterModels::sizeGroup& fi = this->fi();
    const UPtrList<diameterModels::sizeGroup>& popBalFis =
        fi.group().popBal().sizeGroups();

    Pair<tmp<DimensionedField<scalar, volMesh>>> fiSourceCoeffs =
        refCast<const growthFvScalarFieldSource>
        (
            fi.sources()[model.name()]
        ).sourceCoeffs(model);

    Pair<tmp<DimensionedField<scalar, volMesh>>> tsourceCoeffs;

    if (fi.i() != 0)
    {
        tsourceCoeffs.first() = fiSourceCoeffs.first()*value(-1, model);
    }

    if (fi.i() != popBalFis.size() - 1)
    {
        tsourceCoeffs.second() = fiSourceCoeffs.second()*value(+1, model);
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
    const diameterModels::sizeGroup& fi = this->fi();
    const UPtrList<diameterModels::sizeGroup>& velGrpFis =
        fi.group().sizeGroups();

    Pair<tmp<DimensionedField<scalar, volMesh>>> tsourceCoeffs =
        sourceCoeffs(model);

    return
        fi.i() == velGrpFis.first().i()
      ? neg(source)*tsourceCoeffs.second()
      : fi.i() == velGrpFis.last().i()
      ? pos(source)*tsourceCoeffs.first()
      : pos(source)*tsourceCoeffs.first()
      + neg(source)*tsourceCoeffs.second();
}


// ************************************************************************* //

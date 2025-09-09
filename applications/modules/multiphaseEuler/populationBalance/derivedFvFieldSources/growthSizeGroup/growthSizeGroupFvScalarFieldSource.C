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

#include "growthSizeGroupFvScalarFieldSource.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthSizeGroupFvScalarFieldSource::w
(
    const fvSource& model,
    const diameterModels::sizeGroup& fi
) const
{
    // Get the moment
    const label q = this->q();

    // Name of the weight normalisation field
    const word wName = fi.group().phase().name() + ":" + model.name() + ":w";

    // Quick return for volume moments that do not need to compute a weight
    if (q == 3)
    {
        return DimensionedField<scalar, volMesh>::New(wName, fi.mesh(), fi.x());
    }

    // Create the weight normalisation field if it does not yet exist
    const bool haveW = db().foundObject<volScalarField::Internal>(wName);
    if (!haveW)
    {
        volScalarField::Internal* wPtr =
            new volScalarField::Internal
            (
                IOobject
                (
                    wName,
                    internalField().mesh().time().name(),
                    internalField().mesh()
                ),
                internalField().mesh(),
                pow(dimVolume, scalar(q)/3 - 1)
            );

        wPtr->store();
    }

    // Update the weight normalisation field if it is out of date
    volScalarField::Internal& w =
        db().lookupObjectRef<volScalarField::Internal>(wName);
    if (!haveW || !w.hasStoredOldTimes())
    {
        const UPtrList<diameterModels::sizeGroup>& velGrpFis =
            fi.group().sizeGroups();

        w.primitiveFieldRef() = scalar(0);

        forAll(velGrpFis, i)
        {
            w.primitiveFieldRef() +=
                velGrpFis[i]*pow(velGrpFis[i].x(), scalar(q)/3 - 1);
        }
    }

    // Return the normalised weight for this size-group
    return pow(fi.x(), scalar(q)/3)/w;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthSizeGroupFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const diameterModels::sizeGroup& fi =
        refCast<const diameterModels::sizeGroup>(this->internalField());
    const UPtrList<diameterModels::sizeGroup>& popBalFis =
        fi.group().popBal().sizeGroups();

    const dimensionedScalar& xi = fi.x();
    const DimensionedField<scalar, volMesh> wi(w(model, fi));

    tmp<DimensionedField<scalar, volMesh>> tinternalCoeff;

    if (fi.i() == 0)
    {
        tinternalCoeff = neg(source)*wi/xi;
    }
    else
    {
        const dimensionedScalar& xiMinus1 = popBalFis[fi.i() - 1].x();
        tinternalCoeff = neg(source)*wi/(xi - xiMinus1);
    }

    if (fi.i() != popBalFis.size() - 1)
    {
        const dimensionedScalar& xiPlus1 = popBalFis[fi.i() + 1].x();
        tinternalCoeff.ref() -= pos(source)*wi/(xiPlus1 - xi);
    }
    else
    {
        tinternalCoeff.ref() += pos(source)*wi/xi;
    }

    return tinternalCoeff;
}


Foam::Pair<Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>>
Foam::growthSizeGroupFvScalarFieldSource::sourceCoeffs
(
    const fvSource& model
) const
{
    const diameterModels::sizeGroup& fi =
        refCast<const diameterModels::sizeGroup>(this->internalField());
    const UPtrList<diameterModels::sizeGroup>& popBalFis =
        fi.group().popBal().sizeGroups();

    const dimensionedScalar& xi = fi.x();

    Pair<tmp<DimensionedField<scalar, volMesh>>> tsourceCoeffs;

    if (fi.i() != 0)
    {
        const diameterModels::sizeGroup& fiMinus1 = popBalFis[fi.i() - 1];
        const dimensionedScalar& xiMinus1 = fiMinus1.x();
        tsourceCoeffs.first() =
            fiMinus1()*w(model, fiMinus1)*(xi/xiMinus1)/(xi - xiMinus1);
    }

    if (fi.i() != popBalFis.size() - 1)
    {
        const diameterModels::sizeGroup& fiPlus1 = popBalFis[fi.i() + 1];
        const dimensionedScalar& xiPlus1 = fiPlus1.x();
        tsourceCoeffs.second() =
            - fiPlus1()*w(model, fiPlus1)*(xi/xiPlus1)/(xiPlus1 - xi);
    }

    return tsourceCoeffs;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthSizeGroupFvScalarFieldSource::sourceCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const diameterModels::sizeGroup& fi =
        refCast<const diameterModels::sizeGroup>(this->internalField());
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

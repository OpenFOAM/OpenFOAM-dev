/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "growthGroupFractionFvScalarFieldSource.H"
#include "populationBalanceModel.H"
#include "sizeSpecificMassTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        growthGroupFractionFvScalarFieldSource
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::growthGroupFractionFvScalarFieldSource::q
(
    const fvSource& model
) const
{
    return -labelMax;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::growthGroupFractionFvScalarFieldSource::w
(
    const fvSource& model,
    const label i
) const
{
    const populationBalanceModel& popBal = this->popBal();

    // Get the moment
    const label q = this->q(model);

    // Name of the weight normalisation field
    const word wName = popBal.phases()[i].name() + ":" + model.name() + ":w";

    // Quick return for volume moments that do not need to compute a weight
    if (q == 3)
    {
        return
            DimensionedField<scalar, fvMesh>::New
            (
                wName,
                internalField().mesh(),
                popBal.v(i)
            );
    }

    // Create the weight normalisation field if it does not yet exist
    const bool haveW = db().foundObject<volInternalScalarField>(wName);
    if (!haveW)
    {
        volInternalScalarField* wPtr =
            new volInternalScalarField
            (
                IOobject
                (
                    wName,
                    internalField().mesh().time().name(),
                    internalField().mesh()
                ),
                internalField().mesh(),
                pow(dimensions::volume, scalar(q)/3 - 1),
                false
            );

        wPtr->store();
    }

    // Update the weight normalisation field if it is out of date
    volInternalScalarField& w =
        db().lookupObjectRef<volInternalScalarField>(wName);
    if (!haveW || !w.hasStoredOldTimes())
    {
        w.primitiveFieldRef() = scalar(0);

        for
        (
            label j = popBal.diameters()[i].iFirst();
            j <= popBal.diameters()[i].iLast();
            ++ j
        )
        {
            w.primitiveFieldRef() +=
                popBal.f(j).primitiveField()
               *pow(popBal.v(j).value(), scalar(q)/3 - 1);
        }
    }

    // Return the normalised weight for this group
    return pow(popBal.v(i), scalar(q)/3)/w/popBal.v(i);
}


void Foam::growthGroupFractionFvScalarFieldSource::check
(
    const fvSource& model
) const
{
    const bool isSizeSpecific = isA<const fv::sizeSpecificMassTransfer>(model);

    const label q = this->q(model);

    if (isSizeSpecific == (q != -labelMax))
    {
        FatalErrorInFunction
            << "Condition of type " << type() << " cannot be used for source "
            << model.name() << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << " as "
            << (isSizeSpecific ? "both" : "neither") << " the model "
            << (isSizeSpecific ? "and" : "nor") << " the condition define how "
            << "the source is distributed across the groups"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::growthGroupFractionFvScalarFieldSource::
growthGroupFractionFvScalarFieldSource
(
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    growthFvScalarFieldSource(iF, dict),
    groupPropertyFvScalarField(iF)
{}


Foam::growthGroupFractionFvScalarFieldSource::
growthGroupFractionFvScalarFieldSource
(
    const growthGroupFractionFvScalarFieldSource& field,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    growthFvScalarFieldSource(field, iF),
    groupPropertyFvScalarField(iF)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::growthGroupFractionFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, fvMesh>& source
) const
{
    check(model);

    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();

    const dimensionedScalar& xi = popBal.v(i);

    tmp<DimensionedField<scalar, fvMesh>> tinternalCoeff;

    if (i == 0)
    {
        tinternalCoeff = neg(source);
    }
    else
    {
        const dimensionedScalar& xiMinus1 = popBal.v(i - 1);
        tinternalCoeff = neg(source)*xi/(xi - xiMinus1);
    }

    if (i != popBal.nGroups() - 1)
    {
        const dimensionedScalar& xiPlus1 = popBal.v(i + 1);
        tinternalCoeff.ref() -= pos(source)*xi/(xiPlus1 - xi);
    }
    else
    {
        tinternalCoeff.ref() += pos(source);
    }

    if (!isA<const fv::sizeSpecificMassTransfer>(model))
    {
        tinternalCoeff.ref() *= w(model, i);
    }

    return tinternalCoeff;
}


Foam::Pair<Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>>
Foam::growthGroupFractionFvScalarFieldSource::sourceCoeffs
(
    const fvSource& model
) const
{
    check(model);

    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();

    const dimensionedScalar& xi = popBal.v(i);

    Pair<tmp<DimensionedField<scalar, fvMesh>>> tsourceCoeffs;

    if (i != 0)
    {
        const DimensionedField<scalar, fvMesh>& fiMinus1 = popBal.f(i - 1);
        const dimensionedScalar& xiMinus1 = popBal.v(i - 1);
        tsourceCoeffs.first() = fiMinus1*xi/(xi - xiMinus1);

        if (!isA<const fv::sizeSpecificMassTransfer>(model))
        {
            tsourceCoeffs.first().ref() *= w(model, i - 1);
        }
    }

    if (i != popBal.nGroups() - 1)
    {
        const DimensionedField<scalar, fvMesh>& fiPlus1 = popBal.f(i + 1);
        const dimensionedScalar& xiPlus1 = popBal.v(i + 1);
        tsourceCoeffs.second() = -fiPlus1*xi/(xiPlus1 - xi);

        if (!isA<const fv::sizeSpecificMassTransfer>(model))
        {
            tsourceCoeffs.second().ref() *= w(model, i + 1);
        }
    }

    return tsourceCoeffs;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::growthGroupFractionFvScalarFieldSource::sourceTerm
(
    const fvSource& model,
    const DimensionedField<scalar, fvMesh>& source
) const
{
    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();

    Pair<tmp<DimensionedField<scalar, fvMesh>>> tsourceCoeffs =
        sourceCoeffs(model);

    if (isA<const fv::sizeSpecificMassTransfer>(model))
    {
        const fv::sizeSpecificMassTransfer& ssmtModel =
            refCast<const fv::sizeSpecificMassTransfer>(model);

        return
            i == popBal.diameters()[i].iFirst()
          ? eval(negPart(ssmtModel.groupMDotByF(i + 1))*tsourceCoeffs.second())
          : i == popBal.diameters()[i].iLast()
          ? eval(posPart(ssmtModel.groupMDotByF(i - 1))*tsourceCoeffs.first())
          : eval
            (
                posPart(ssmtModel.groupMDotByF(i - 1))*tsourceCoeffs.first()
              + negPart(ssmtModel.groupMDotByF(i + 1))*tsourceCoeffs.second()
            );
    }
    else
    {
        return
            i == popBal.diameters()[i].iFirst()
          ? eval(negPart(source)*tsourceCoeffs.second())
          : i == popBal.diameters()[i].iLast()
          ? eval(posPart(source)*tsourceCoeffs.first())
          : eval
            (
                posPart(source)*tsourceCoeffs.first()
              + negPart(source)*tsourceCoeffs.second()
            );
    }
}


// ************************************************************************* //

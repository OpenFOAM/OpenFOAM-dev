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

#include "growthGroupFractionFvScalarFieldSource.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthGroupFractionFvScalarFieldSource::w
(
    const fvSource& model,
    const label i
) const
{
    const populationBalanceModel& popBal = this->popBal();

    // Get the moment
    const label q = this->q();

    // Name of the weight normalisation field
    const word wName = popBal.phases()[i].name() + ":" + model.name() + ":w";

    // Quick return for volume moments that do not need to compute a weight
    if (q == 3)
    {
        return
            DimensionedField<scalar, volMesh>::New
            (
                wName,
                internalField().mesh(),
                popBal.v(i)
            );
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
        w.primitiveFieldRef() = scalar(0);

        for
        (
            label j = popBal.diameters()[i].iFirst();
            j <= popBal.diameters()[i].iLast();
            ++ j
        )
        {
            w.primitiveFieldRef() +=
                popBal.f(j)*pow(popBal.v(j), scalar(q)/3 - 1);
        }
    }

    // Return the normalised weight for this group
    return pow(popBal.v(i), scalar(q)/3)/w;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::growthGroupFractionFvScalarFieldSource::
growthGroupFractionFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
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
    const DimensionedField<scalar, volMesh>& iF
)
:
    growthFvScalarFieldSource(field, iF),
    groupPropertyFvScalarField(iF)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthGroupFractionFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();

    const dimensionedScalar& xi = popBal.v(i);
    const DimensionedField<scalar, volMesh> wi(w(model, i));

    tmp<DimensionedField<scalar, volMesh>> tinternalCoeff;

    if (i == 0)
    {
        tinternalCoeff = neg(source)*wi/xi;
    }
    else
    {
        const dimensionedScalar& xiMinus1 = popBal.v(i - 1);
        tinternalCoeff = neg(source)*wi/(xi - xiMinus1);
    }

    if (i != popBal.nGroups() - 1)
    {
        const dimensionedScalar& xiPlus1 = popBal.v(i + 1);
        tinternalCoeff.ref() -= pos(source)*wi/(xiPlus1 - xi);
    }
    else
    {
        tinternalCoeff.ref() += pos(source)*wi/xi;
    }

    return tinternalCoeff;
}


Foam::Pair<Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>>
Foam::growthGroupFractionFvScalarFieldSource::sourceCoeffs
(
    const fvSource& model
) const
{
    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();

    const dimensionedScalar& xi = popBal.v(i);

    Pair<tmp<DimensionedField<scalar, volMesh>>> tsourceCoeffs;

    if (i != 0)
    {
        const DimensionedField<scalar, volMesh>& fiMinus1 = popBal.f(i - 1);
        const dimensionedScalar& xiMinus1 = popBal.v(i - 1);
        tsourceCoeffs.first() =
            fiMinus1*w(model, i - 1)*(xi/xiMinus1)/(xi - xiMinus1);
    }

    if (i != popBal.nGroups() - 1)
    {
        const DimensionedField<scalar, volMesh>& fiPlus1 = popBal.f(i + 1);
        const dimensionedScalar& xiPlus1 = popBal.v(i + 1);
        tsourceCoeffs.second() =
            -fiPlus1*w(model, i + 1)*(xi/xiPlus1)/(xiPlus1 - xi);
    }

    return tsourceCoeffs;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthGroupFractionFvScalarFieldSource::sourceCoeff
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
        i == popBal.diameters()[i].iFirst()
      ? neg(source)*tsourceCoeffs.second()
      : i == popBal.diameters()[i].iLast()
      ? pos(source)*tsourceCoeffs.first()
      : pos(source)*tsourceCoeffs.first()
      + neg(source)*tsourceCoeffs.second();
}


// ************************************************************************* //

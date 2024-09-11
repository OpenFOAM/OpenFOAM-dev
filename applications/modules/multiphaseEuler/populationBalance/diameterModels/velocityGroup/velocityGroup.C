/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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

#include "velocityGroup.H"
#include "addToRunTimeSelectionTable.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(velocityGroup, 0);
    addToRunTimeSelectionTable(diameterModel, velocityGroup, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::dsm() const
{
    tmp<volScalarField> tsumFi
    (
        volScalarField::New
        (
            "sumFi",
            phase().mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );
    tmp<volScalarField> tsumFiAbyV
    (
        volScalarField::New
        (
            "sumFiAbyV",
            phase().mesh(),
            dimensionedScalar(inv(dimLength), Zero)
        )
    );

    volScalarField& sumFi = tsumFi.ref();
    volScalarField& sumFiAbyV = tsumFiAbyV.ref();

    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = sizeGroups_[i];

        sumFi += max(fi, rootVSmall);
        sumFiAbyV += max(fi, rootVSmall)*fi.a()/fi.x();
    }

    return 6*sumFi/tsumFiAbyV;
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::fSum() const
{
    tmp<volScalarField> tsumSizeGroups
    (
        volScalarField::New
        (
            "sumSizeGroups",
            phase().mesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField& sumSizeGroups = tsumSizeGroups.ref();

    forAll(sizeGroups_, i)
    {
        sumSizeGroups += sizeGroups_[i];
    }

    return tsumSizeGroups;
}


void Foam::diameterModels::velocityGroup::scale()
{
    Info<< "Scaling sizeGroups for velocityGroup " << phase().name() << endl;

    forAll(sizeGroups_, i)
    {
        sizeGroups_[i].max(0);
    };

    const volScalarField fSum(this->fSum());

    forAll(sizeGroups_, i)
    {
        sizeGroups_[i] /= fSum;

        sizeGroups_[i].correctBoundaryConditions();
    };
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::velocityGroup::velocityGroup
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    popBalName_(diameterProperties.lookup("populationBalance")),
    popBalPtr_(nullptr),
    sizeGroups_
    (
        diameterProperties.lookup("sizeGroups"),
        sizeGroup::iNew
        (
            *this,
            populationBalanceModel::groups::New
            (
                popBalName_,
                phase.mesh()
            ).nSizeGroups()
        )
    ),
    d_(IOobject::groupName("d", phase.name()), dsm())
{
    populationBalanceModel::groups::New
    (
        popBalName_,
        phase.mesh()
    ).insert(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::velocityGroup::~velocityGroup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::diameterModels::populationBalanceModel&
Foam::diameterModels::velocityGroup::popBal() const
{
    if (popBalPtr_ == nullptr)
    {
        popBalPtr_ =
            &phase().mesh().lookupObject<populationBalanceModel>(popBalName_);
    }

    return *popBalPtr_;
}


Foam::tmp<Foam::volScalarField> Foam::diameterModels::velocityGroup::d() const
{
    return d_;
}


Foam::tmp<Foam::volScalarField> Foam::diameterModels::velocityGroup::Av() const
{
    tmp<volScalarField> tsumFiAbyV
    (
        volScalarField::New
        (
            "sumFiAbyV",
            phase().mesh(),
            dimensionedScalar(inv(dimLength), Zero)
        )
    );
    volScalarField& sumFiAbyV = tsumFiAbyV.ref();

    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = sizeGroups_[i];

        sumFiAbyV += fi*fi.a()/fi.x();
    }

    return phase()*tsumFiAbyV;
}


void Foam::diameterModels::velocityGroup::correct()
{
    const populationBalanceModel& popBal = this->popBal();

    if (!popBal.solveOnFinalIterOnly() || popBal.fluid().pimple().finalIter())
    {
        forAll(sizeGroups_, i)
        {
            sizeGroups_[i].correct();
        }

        if
        (
            phase()
           .mesh()
           .solution()
           .solverDict(popBalName_)
           .lookupOrDefault<Switch>
            (
                "scale",
                true
            )
        )
        {
            scale();
        }

        volScalarField::Internal fSum(this->fSum());

        Info<< phase().name() << " sizeGroups-sum volume fraction, min, max = "
            << fSum.weightedAverage(phase().mesh().V()).value()
            << ' ' << min(fSum).value()
            << ' ' << max(fSum).value()
            << endl;

        d_ = dsm();

        Info<< this->phase().name() << " Sauter mean diameter, min, max = "
            << d_.weightedAverage(d_.mesh().V()).value()
            << ' ' << min(d_).value()
            << ' ' << max(d_).value()
            << endl;
    }
}


bool Foam::diameterModels::velocityGroup::read
(
    const dictionary& phaseProperties
)
{
    diameterModel::read(phaseProperties);

    return true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
#include "sizeGroup.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

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

Foam::tmp<Foam::volScalarField> Foam::diameterModels::velocityGroup::dsm() const
{
    tmp<volScalarField> tInvDsm
    (
        volScalarField::New
        (
            "invDsm",
            phase().mesh(),
            dimensionedScalar(inv(dimLength), Zero)
        )
    );

    volScalarField& invDsm = tInvDsm.ref();

    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = sizeGroups_[i];

        invDsm += fi.a()*fi/fi.x();
    }

    return 6.0/tInvDsm;
}


Foam::tmp<Foam::volScalarField> Foam::diameterModels::velocityGroup::N() const
{
    tmp<volScalarField> tN
    (
        volScalarField::New
        (
            "N",
            phase().mesh(),
            dimensionedScalar(inv(dimVolume), 0)
        )
    );

    volScalarField& N = tN.ref();

    forAll(sizeGroups_, i)
    {
        N += phase()*sizeGroups_[i]/sizeGroups_[i].x();
    }

    return tN;
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

    f_ = fSum();

    forAll(sizeGroups_, i)
    {
        sizeGroups_[i] /= f_;

        sizeGroups_[i].correctBoundaryConditions();
    };
}


// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::calcD() const
{
    return d_;
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::calcA() const
{
    tmp<volScalarField> tA
    (
        volScalarField::New
        (
            "a",
            phase().mesh(),
            dimensionedScalar(inv(dimLength), Zero)
        )
    );

    volScalarField& a = tA.ref();

    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = sizeGroups_[i];

        a += fi.a()*fi/fi.x();
    }

    return phase()*a;
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
    f_
    (
        IOobject
        (
            IOobject::groupName
            (
                "f",
                IOobject::groupName
                (
                    phase.name(),
                    popBalName_
                )
            ),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    sizeGroups_
    (
        diameterProperties.lookup("sizeGroups"),
        sizeGroup::iNew(phase, *this)
    ),
    d_(dRef())
{
    d_ = dsm();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::velocityGroup::~velocityGroup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::velocityGroup::correct()
{
    forAll(sizeGroups_, i)
    {
        sizeGroups_[i].correct();
    }

    if
    (
        phase().mesh().solverDict(popBalName_).lookupOrDefault<Switch>
        (
            "scale",
            true
        )
    )
    {
        scale();
    }

    f_ = fSum();

    f_.correctBoundaryConditions();

    Info<< phase().name() << " sizeGroups-sum volume fraction, min, max = "
        << f_.weightedAverage(phase().mesh().V()).value()
        << ' ' << min(f_).value()
        << ' ' << max(f_).value()
        << endl;

    d_ = dsm();

    Info<< this->phase().name() << " Sauter mean diameter, min, max = "
        << d_.weightedAverage(d_.mesh().V()).value()
        << ' ' << min(d_).value()
        << ' ' << max(d_).value()
        << endl;
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

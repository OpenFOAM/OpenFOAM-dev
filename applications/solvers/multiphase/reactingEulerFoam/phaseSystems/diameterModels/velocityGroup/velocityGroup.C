/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

    addToRunTimeSelectionTable
    (
        diameterModel,
        velocityGroup,
        dictionary
    );
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
            phase_.mesh(),
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
            phase_.mesh(),
            dimensionedScalar(inv(dimVolume), 0)
        )
    );

    volScalarField& N = tN.ref();

    forAll(sizeGroups_, i)
    {
        N += phase_*sizeGroups_[i]/sizeGroups_[i].x();
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
            phase_.mesh(),
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
    Info<< "Scaling sizeGroups for velocityGroup " << phase_.name() << endl;

    forAll(sizeGroups_, i)
    {
        sizeGroups_[i].max(0);
    };

    f_ = fSum();

    forAll(sizeGroups_, i)
    {
        sizeGroups_[i] /= f_;
    };
}


Foam::tmp<Foam::fv::convectionScheme<Foam::scalar>>
Foam::diameterModels::velocityGroup::mvconvection() const
{
    tmp<fv::convectionScheme<Foam::scalar>> mvConvection
    (
        fv::convectionScheme<Foam::scalar>::New
        (
            phase_.mesh(),
            fields_,
            phase_.alphaRhoPhi(),
            phase_.mesh().divScheme
            (
                "div(" + phase_.alphaRhoPhi()().name() + ",f)"
            )
        )
    );

    return mvConvection;
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
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    dmdt_
    (
        IOobject
        (
            IOobject::groupName("source", phase.name()),
            phase.time().timeName(),
            phase.mesh()
        ),
        phase.mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    )
{
    forAll(sizeGroups_, i)
    {
        fields_.add(sizeGroups_[i]);
    }

    d_ = dsm();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::velocityGroup::~velocityGroup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::diameterModels::velocityGroup::preSolve()
{
    mvConvection_ = mvconvection();
}


void Foam::diameterModels::velocityGroup::postSolve()
{
    if
    (
        phase_.mesh().solverDict(popBalName_).lookupOrDefault<Switch>
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

    Info<< phase_.name() << " sizeGroups-sum volume fraction, min, max = "
        << f_.weightedAverage(phase_.mesh().V()).value()
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


bool Foam::diameterModels::velocityGroup::
read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::d() const
{
    return d_;
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::secondMoment() const
{
    tmp<volScalarField> tm2
    (
        new volScalarField
        (
            IOobject
            (
                "m2",
                phase_.time().timeName(),
                phase_.mesh()
            ),
            phase_.mesh(),
            dimensionedScalar("m2", inv(dimLength), Zero)
        )
    );

    volScalarField& m2 = tm2.ref();

    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = sizeGroups_[i];

        m2 += sqr(fi.d())*formFactor()*fi
           *max(fi.phase(), small)/fi.x();
    }

    return tm2;
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::thirdMoment() const
{
    tmp<volScalarField> tm3
    (
        new volScalarField
        (
            IOobject
            (
                "m3",
                phase_.time().timeName(),
                phase_.mesh()
            ),
            phase_.mesh(),
            dimensionedScalar("m3", dimless, Zero)
        )
    );

    volScalarField& m3 = tm3.ref();

    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = sizeGroups_[i];

        m3 += pow3(fi.d())*formFactor()*fi
           *max(fi.phase(), small)/fi.x();
    }

    return tm3;
}


Foam::tmp<Foam::volScalarField> Foam::diameterModels::velocityGroup::dsm() const
{
    return
        max(min(phase_/m2_, sizeGroups_.last().d()), sizeGroups_.first().d());
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::velocityGroup::fSum() const
{
    tmp<volScalarField> tsumSizeGroups
    (
        new volScalarField
        (
            IOobject
            (
                "sumSizeGroups",
                phase_.time().timeName(),
                phase_.mesh()
            ),
            phase_.mesh(),
            dimensionedScalar("sumSizeGroups", dimless, 0)
        )
    );

    volScalarField& sumSizeGroups = tsumSizeGroups.ref();

    forAll(sizeGroups_, i)
    {
        sumSizeGroups += sizeGroups_[i];
    }

    return tsumSizeGroups;
}


void Foam::diameterModels::velocityGroup::renormalize()
{
        Info<< phase_.name()
            << " renormalizing sizeGroups"
            << endl;

        // Set negative values to zero
        forAll(sizeGroups_, i)
        {
            sizeGroups_[i] *= pos(sizeGroups_[i]);
        };

        forAll(sizeGroups_, i)
        {
            sizeGroups_[i] /= fSum_;
        };
}


Foam::tmp<Foam::fv::convectionScheme<Foam::scalar>>
Foam::diameterModels::velocityGroup::mvconvection() const
{
    tmp<fv::convectionScheme<Foam::scalar> > mvConvection
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
    formFactor_("formFactor", dimless, diameterProperties),
    sizeGroups_
    (
        diameterProperties.lookup("sizeGroups"),
        sizeGroup::iNew(phase, *this)
    ),
    fSum_
    (
        IOobject
        (
            IOobject::groupName
            (
                "fsum",
                IOobject::groupName
                (
                    phase.name(),
                    popBalName_
                )
            ),
            phase.time().timeName(),
            phase.mesh()
        ),
        fSum()
    ),
    m2_
    (
        IOobject
        (
            IOobject::groupName
            (
                "m2",
                IOobject::groupName
                (
                    phase.name(),
                    popBalName_
                )
            ),
            phase.time().timeName(),
            phase.mesh()
        ),
        phase.mesh(),
        dimensionedScalar("m2", inv(dimLength), Zero)
    ),
    m3_
    (
        IOobject
        (
            IOobject::groupName
            (
                "m3",
                IOobject::groupName
                (
                    phase.name(),
                    popBalName_
                )
            ),
            phase.time().timeName(),
            phase.mesh()
        ),
        phase.mesh(),
        dimensionedScalar("m3", dimless, Zero)
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
        dimensionedScalar("d", dimLength, Zero)
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
        dimensionedScalar("dmdt", dimDensity/dimTime, Zero)
    )
{
    if
    (
        phase_.mesh().solverDict(popBalName_).lookupOrDefault<Switch>
        (
            "renormalizeAtRestart",
            false
        )
     ||
        phase_.mesh().solverDict(popBalName_).lookupOrDefault<Switch>
        (
            "renormalize",
            false
        )
    )
    {
        renormalize();
    }

    fSum_ = fSum();

    if
    (
        mag(1 - fSum_.weightedAverage(fSum_.mesh().V()).value()) >= 1e-5
     || mag(1 - max(fSum_).value()) >= 1e-5
     || mag(1 - min(fSum_).value()) >= 1e-5
    )
    {
        FatalErrorInFunction
            << " Initial values of the sizeGroups belonging to velocityGroup "
            << this->phase().name()
            << " must add to" << nl << " unity. This condition might be"
            << " violated due to wrong entries in the" << nl
            << " velocityGroupCoeffs subdictionary or bad initial conditions in"
            << " the startTime" << nl
            << " directory. The sizeGroups can be renormalized at every"
            << " timestep or at restart" << nl
            << " only by setting the corresponding switch renormalize or"
            << " renormalizeAtRestart" << nl
            << " in the fvSolution subdictionary " << popBalName_ << "."
            << " Note that boundary conditions are not" << nl << "renormalized."
            << exit(FatalError);
    }

    forAll(sizeGroups_, i)
    {
        fields_.add(sizeGroups_[i]);
    }

    m2_ = secondMoment();

    m3_ = thirdMoment();

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
    m2_ = secondMoment();

    m3_ = thirdMoment();

    d_ = dsm();

    Info<< this->phase().name() << " Sauter mean diameter, min, max = "
        << d_.weightedAverage(d_.mesh().V()).value()
        << ' ' << min(d_).value()
        << ' ' << max(d_).value()
        << endl;

    fSum_ = fSum();

    Info<< phase_.name() << " sizeGroups-sum volume fraction, min, max = "
        << fSum_.weightedAverage(phase_.mesh().V()).value()
        << ' ' << min(fSum_).value()
        << ' ' << max(fSum_).value()
        << endl;

    if
    (
        phase_.mesh().solverDict(popBalName_).lookupOrDefault<Switch>
        (
            "renormalize",
            false
        )
    )
    {
        renormalize();
    }
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

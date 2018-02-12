/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "IATE.H"
#include "IATEsource.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcAverage.H"
#include "fvOptions.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(IATE, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        IATE,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATE::IATE
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    kappai_
    (
        IOobject
        (
            IOobject::groupName("kappai", phase.name()),
            phase_.time().timeName(),
            phase_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.mesh()
    ),
    dMax_("dMax", dimLength, diameterProperties_),
    dMin_("dMin", dimLength, diameterProperties_),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        diameterProperties_
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase_.time().timeName(),
            phase_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dsm()
    ),
    sources_
    (
        diameterProperties_.lookup("sources"),
        IATEsource::iNew(*this)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::IATE::~IATE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::IATE::dsm() const
{
    return max(6/max(kappai_, 6/dMax_), dMin_);
}

void Foam::diameterModels::IATE::correct()
{
    volScalarField alphaAv
    (
        max
        (
            0.5*fvc::average(phase_ + phase_.oldTime()),
            residualAlpha_
        )
    );

    // Initialise the accumulated source term to the dilatation effect
    fvScalarMatrix R
    (
       -fvm::SuSp
        (
            ((1.0/3.0)/alphaAv)
           *(
                fvc::ddt(phase_) + fvc::div(phase_.alphaPhi())
              - phase_.continuityError()/phase_.rho()
            ),
            kappai_
        )
    );

    // Accumulate the run-time selectable sources
    forAll(sources_, j)
    {
        R += sources_[j].R(alphaAv, kappai_);
    }

    fv::options& fvOptions(fv::options::New(phase_.mesh()));

    // Construct the interfacial curvature equation
    fvScalarMatrix kappaiEqn
    (
        fvm::ddt(kappai_) + fvm::div(phase_.phi(), kappai_)
      - fvm::Sp(fvc::div(phase_.phi()), kappai_)
     ==
        R
      + fvOptions(kappai_)
    );

    kappaiEqn.relax();

    fvOptions.constrain(kappaiEqn);

    kappaiEqn.solve();

    // Update the Sauter-mean diameter
    d_ = dsm();
}


bool Foam::diameterModels::IATE::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    diameterProperties_.lookup("dMax") >> dMax_;
    diameterProperties_.lookup("dMin") >> dMin_;

    // Re-create all the sources updating number, type and coefficients
    PtrList<IATEsource>
    (
        diameterProperties_.lookup("sources"),
        IATEsource::iNew(*this)
    ).transfer(sources_);

    return true;
}


// ************************************************************************* //

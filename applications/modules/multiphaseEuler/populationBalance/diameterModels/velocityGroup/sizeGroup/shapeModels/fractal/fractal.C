/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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

#include "fractal.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{
    defineTypeNameAndDebug(fractal, 0);
    addToRunTimeSelectionTable(shapeModel, fractal, dictionary);
}
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::fractal::dColl() const
{
    tmp<volScalarField> tDColl
    (
        volScalarField::New
        (
            "dColl",
            group().mesh(),
            dimensionedScalar(dimLength, Zero)
        )
    );

    volScalarField& dColl = tDColl.ref();

    dColl =
        6/kappa_
       *pow(group().x()*pow3(kappa_)/(36*pi*alphaC_), 1/Df_);

    return tDColl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::fractal::fractal
(
    const dictionary& dict,
    const sizeGroup& group,
    const dictionary& groupDict
)
:
    SecondaryPropertyModel<shapeModel>(group),
    kappa_
    (
        sizeGroup::fieldIo("kappa", group.i(), group.group()),
        sizeGroup::field("kappa", group.i(), group.group())
    ),
    Df_("Df", dimless, groupDict),
    alphaC_("alphaC", dimless, groupDict),
    dColl_
    (
        IOobject
        (
            "dColl" + group.name().substr(1),
            group.mesh().time().name(),
            group.mesh()
        ),
        this->dColl()
    ),
    Su_
    (
        IOobject
        (
            IOobject::groupName("Su", kappa_.name()),
            group.mesh().time().name(),
            group.mesh()
        ),
        group.mesh(),
        dimensionedScalar(kappa_.dimensions()/dimTime, Zero)
    )
{
    // Check and filter for old syntax (remove in due course)
    if (groupDict.found("kappa"))
    {
        FatalErrorInFunction
            << "A 'kappa' entry should not be specified for size-group #"
            << group.i() << " of population balance "
            << group.group().popBalName()
            << ". Instead, the value should be initialised within the field, "
            << kappa_.name() << " (or the default field, "
            << IOobject::groupName("kappaDefault", group.phase().name())
            << ", as appropriate)."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::fractal::~fractal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::diameterModels::shapeModels::fractal::fld() const
{
    return kappa_;
}


Foam::volScalarField::Internal&
Foam::diameterModels::shapeModels::fractal::src()
{
    return Su_;
}


void Foam::diameterModels::shapeModels::fractal::correct()
{
    const sizeGroup& fi = group();
    const label i = fi.i();
    const phaseModel& phase = fi.phase();
    const volScalarField& alpha = phase;
    const volScalarField& rho = phase.rho();

    const populationBalanceModel& popBal = fi.group().popBal();

    const volScalarField alphaFi
    (
        IOobject::groupName
        (
            alpha.member() + fi.member().capitalise(),
            alpha.group()
        ),
        alpha*fi
    );

    const surfaceScalarField alphaFiPhi
    (
        IOobject::groupName
        (
            alpha.member() + fi.member().capitalise() + "Phi",
            alpha.group()
        ),
        max(fvc::interpolate(fi, "fi"), small)*phase.alphaPhi()
    );

    // !!! Create pointers to other kappa fields. The shortcut taken here is
    // that only adjacent fields are set. We know that the transfers never
    // extend beyond these adjacent groups. Eventually the plan is to store all
    // the fields in a pointer list anyway, rather than having nested per-group
    // models, so then it will be possible just to pass the list directly.
    UPtrList<const volScalarField> flds(popBal.sizeGroups().size());
    for (label deltai = -1; deltai <= +1; ++ deltai)
    {
        const label j = fi.i() + deltai;
        if (j < 0 || j >= popBal.sizeGroups().size()) continue;
        const sizeGroup& fj = popBal.sizeGroups()[j];
        flds.set(fj.i(), &model(fj).fld());
    }

    fvScalarMatrix kappaEqn
    (
        fvm::ddt(alpha, fi, kappa_)
      + fvm::div(alphaFiPhi, kappa_)
     ==
        Su_ + fvm::Sp(popBal.Sp(i)*fi, kappa_)
      + popBal.expansionSu(i, flds) + fvm::Sp(popBal.expansionSp(i)*fi, kappa_)
      + popBal.modelSourceSu(i, flds)
      + popBal.fluid().fvModels().source(alphaFi, rho, kappa_)/rho
      - correction
        (
            fvm::Sp
            (
                max(phase.residualAlpha() - alpha*fi, scalar(0))
               /kappa_.mesh().time().deltaT(),
                kappa_
            )
        )
    );

    kappaEqn.relax();

    popBal.fluid().fvConstraints().constrain(kappaEqn);

    kappaEqn.solve();

    popBal.fluid().fvConstraints().constrain(kappa_);

    // Bound kappa so that the surface-area-volume ratio is greater than that
    // of spherical particles of this group, but less than that of the
    // particles represented by the first size group
    const sizeGroup& f0 = popBal.sizeGroups().first();
    kappa_ = min(max(kappa_, 6/fi.dSph()), 6/f0.dSph());

    kappa_.correctBoundaryConditions();

    // Update the collisional diameter
    dColl_ = dColl();
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::fractal::a() const
{
    return kappa_*group().x();
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::fractal::d() const
{
    return dColl_;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "sinteringModel.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "mixedFvPatchField.H"

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

template<>
const char*
Foam::NamedEnum
<
    Foam::diameterModels::shapeModels::fractal::surfaceGrowthTypes,
    4
>::names[] = {"unknown", "hardSphere", "ParkRogak", "conserved"};

const Foam::NamedEnum
<
    Foam::diameterModels::shapeModels::fractal::surfaceGrowthTypes,
    4
> Foam::diameterModels::shapeModels::fractal::sgTypeNames_;


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
    ),
    sinteringModel_
    (
        sinteringModel::New(dict.subDict(type() + "Coeffs"), *this)
    ),
    sgType_
    (
        sgTypeNames_
        [
            dict.subDict(type() + "Coeffs")
           .lookupOrDefault<word>("surfaceGrowthType", sgTypeNames_.names[0])
        ]
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
            << this->name() << " (or the default field, "
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


Foam::volScalarField&
Foam::diameterModels::shapeModels::fractal::src()
{
    return Su_;
}


void Foam::diameterModels::shapeModels::fractal::correct()
{
    const sizeGroup& fi = group();
    const phaseModel& phase = fi.phase();
    const volScalarField& alpha = phase;

    const populationBalanceModel& popBal =
        group().mesh().lookupObject<populationBalanceModel>
        (
            group().group().popBalName()
        );

    surfaceScalarField fAlphaPhi
    (
        "fAlphaPhi",
        max(fvc::interpolate(fi, "fi"), small)*phase.alphaPhi()
    );

    fvScalarMatrix kappaEqn
    (
        fvm::ddt(alpha, fi, kappa_)
      + fvm::div(fAlphaPhi, kappa_)
      ==
      - sinteringModel_->R()
      + Su_
      - fvm::Sp(popBal.Sp(fi.i())*fi, kappa_)
      - correction
        (
            fvm::Sp
            (
                max(phase.residualAlpha() - alpha*fi, scalar(0))
               /group().mesh().time().deltaT(),
                kappa_
            )
        )
    );

    kappaEqn.relax();

    kappaEqn.solve();

    // Bounding of kappa assuming first sizeGroup to represent one primary
    // particle
    kappa_ =
        min
        (
            max(kappa_, 6/group().dSph()),
            6/popBal.sizeGroups().first().dSph()
        );

    kappa_.correctBoundaryConditions();

    // Update the collisional diameter
    dColl_ = dColl();
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::fractal::a() const
{
    return kappa_*group().x();
}


void Foam::diameterModels::shapeModels::fractal::addDrift
(
    const volScalarField &Su,
    const sizeGroup &fu,
    const driftModel &model
)
{
    const volScalarField& sourceKappa =
        SecondaryPropertyModelTable()[SecondaryPropertyName(fu)]->fld();

    switch (sgType_)
    {
        case sgHardSphere:
        {
            Su_ += sourceKappa*fu.dSph()/group().dSph()*Su;

            break;
        }

        case sgParkRogak:
        {
            const fractal& sourceShape =
                refCast<const fractal>
                (
                    fu.shapeModelPtr()()
                );

            volScalarField dp(6/sourceKappa);
            const volScalarField a(sourceKappa*fu.x());
            const dimensionedScalar dv(group().x() - fu.x());

            const volScalarField da1
            (
                (2.0/3.0)*dv
               *(
                    sourceKappa
                  + sourceShape.Df_*(1/sourceShape.d() - 1/dp)
                )
            );

            dp += 6*(dv*a - fu.x()*da1)/sqr(a);

            const volScalarField np(6*group().x()/pi/pow3(dp));
            const volScalarField dc(dp*pow(np/alphaC_, 1/Df_));

            const volScalarField da2
            (
                dv*(4/dp + 2*Df_/3*(1/dc - 1/dp))
            );

            Su_ += (a + 0.5*da1 + 0.5*da2)/group().x()*Su;

            break;
        }

        case sgConserved:
        {
            SecondaryPropertyModel::addDrift(Su, fu, model);

            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unknown surface growth type. Valid types are:"
                << sgTypeNames_ << nl << exit(FatalError);
        }
    }
}


// ************************************************************************* //

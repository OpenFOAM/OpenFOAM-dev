/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
#include "sinteringModel.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvcSup.H"
#include "fvcDdt.H"
#include "mixedFvPatchField.H"
#include "mathematicalConstants.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

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
    3
>::names[] = {"hardSphere", "ParkRogak", "conserved"};

const Foam::NamedEnum
<
    Foam::diameterModels::shapeModels::fractal::surfaceGrowthTypes,
    3
> Foam::diameterModels::shapeModels::fractal::sgTypeNames_;


using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::fractal::fractal
(
    const dictionary& dict,
    const sizeGroup& group
)
:
    SecondaryPropertyModel<shapeModel>(dict, group),
    kappa_
    (
        IOobject
        (
            "kappa" + group.name().substr(1),
            group.mesh().time().timeName(),
            group.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        group.mesh(),
        dimensionedScalar
        (
            "kappa",
            inv(dimLength),
            group.dict()
        ),
        group.VelocityGroup().f().boundaryField().types()
    ),
    Df_("Df", dimless, group.dict()),
    alphaC_("alphaC", dimless, group.dict()),
    dColl_
    (
        IOobject
        (
            "dColl" + group.name().substr(1),
            group.mesh().time().timeName(),
            group.mesh()
        ),
        this->dColl()
    ),
    Su_
    (
        IOobject
        (
            IOobject::groupName("Su", kappa_.name()),
            group.mesh().time().timeName(),
            group.mesh()
        ),
        group.mesh(),
        dimensionedScalar(kappa_.dimensions()/dimTime, Zero)
    )
{
    // Adjust refValue at mixedFvPatchField boundaries
    forAll(kappa_.boundaryField(), patchi)
    {
        typedef mixedFvPatchField<scalar> mixedFvPatchScalarField;

        if
        (
            isA<mixedFvPatchScalarField>(kappa_.boundaryFieldRef()[patchi])
        )
        {
            mixedFvPatchScalarField& kappa =
                refCast<mixedFvPatchScalarField>
                (
                    kappa_.boundaryFieldRef()[patchi]
                );

            kappa.refValue() = sizeGroup_.dict().lookup<scalar>("kappa");
        }
    }

    sinteringModel_ =
        sinteringModel::New(dict.subDict(type() + "Coeffs"), *this);
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


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::fractal::dColl() const
{
    tmp<volScalarField> tDColl
    (
        volScalarField::New
        (
            "dColl",
            sizeGroup_.mesh(),
            dimensionedScalar(dimLength, Zero)
        )
    );

    volScalarField& dColl = tDColl.ref();

    dColl =
        6.0/kappa_
       *pow(sizeGroup_.x()*pow3(kappa_)/(36.0*pi*alphaC_), 1.0/Df_);

    return tDColl;
}


void Foam::diameterModels::shapeModels::fractal::correct()
{
    const sizeGroup& fi = sizeGroup_;
    const phaseModel& phase = fi.phase();
    const volScalarField& alpha = phase;
    const volScalarField& rho = phase.thermo().rho();

    const populationBalanceModel& popBal =
        sizeGroup_.mesh().lookupObject<populationBalanceModel>
        (
            sizeGroup_.VelocityGroup().popBalName()
        );

    surfaceScalarField fAlphaRhoPhi
    (
        "fAlphaRhoPhi",
        max(fvc::interpolate(fi, "fi"), SMALL)*phase.alphaRhoPhi()
    );

    fvScalarMatrix kappaEqn
    (
        fvc::ddt(alpha, rho, fi)*kappa_.oldTime()
      + alpha*rho*fi*fvm::ddt(kappa_)
      + fvm::div(fAlphaRhoPhi, kappa_)
      + fvm::SuSp
        (
            fi
           *(
                fi.VelocityGroup().dmdt()
              - (fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi()))
            ),
            kappa_
        )
      ==
      - sinteringModel_->R()
      + fvc::Su(Su_*rho, kappa_)
      - fvm::SuSp(popBal.SuSp(fi.i()())*fi*rho, kappa_)
      + fvc::ddt(fi.phase().residualAlpha()*rho, kappa_)
      - fvm::ddt(fi.phase().residualAlpha()*rho, kappa_)
    );

    kappaEqn.relax();

    kappaEqn.solve();

    // Bounding of kappa assuming first sizeGroup to represent one primary
    // particle
    kappa_ =
        min
        (
            max(kappa_, 6.0/sizeGroup_.dSph()),
            6.0/popBal.sizeGroups().first().dSph()
        );

    kappa_.correctBoundaryConditions();

    // Update the collisional diameter
    dColl_ = dColl();
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::shapeModels::fractal::a() const
{
    return kappa_*sizeGroup_.x();
}


void Foam::diameterModels::shapeModels::fractal::addDrift
(
    const volScalarField &Su,
    const sizeGroup &fu,
    const driftModel &model
)
{
    surfaceGrowthTypes sgType
    (
        sgTypeNames_.read
        (
            model.dict().lookup("surfaceGrowthType")
        )
    );

    const volScalarField& sourceKappa =
        SecondaryPropertyModelTable()[SecondaryPropertyName(fu)]->fld();

    switch (sgType)
    {
        case sgHardSphere:
        {
            Su_ += sourceKappa*fu.dSph()/sizeGroup_.dSph()*Su;

            break;
        }

        case sgParkRogak:
        {
            const fractal& sourceShape =
                refCast<const fractal>
                (
                    fu.shapeModelPtr()()
                );

            volScalarField dp(6.0/sourceKappa);
            const volScalarField a(sourceKappa*fu.x());
            const dimensionedScalar dv(sizeGroup_.x() - fu.x());

            const volScalarField da1
            (
                (2.0/3.0)*dv
               *(
                    sourceKappa
                  + sourceShape.Df_*(1.0/sourceShape.d() - 1.0/dp)
                )
            );

            dp += 6.0*(dv*a - fu.x()*da1)/sqr(a);

            const volScalarField np(6.0*sizeGroup_.x()/pi/pow3(dp));
            const volScalarField dc(dp*pow(np/alphaC_, 1.0/Df_));

            const volScalarField da2
            (
                dv*(4.0/dp + 2.0*Df_/3.0*(1.0/dc - 1.0/dp))
            );

            Su_ += (a + 0.5*da1 + 0.5*da2)/sizeGroup_.x()*Su;

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

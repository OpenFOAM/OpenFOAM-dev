/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "FSD.H"
#include "addToRunTimeSelectionTable.H"
#include "LESModel.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(FSD, 0);
    addToRunTimeSelectionTable(combustionModel, FSD, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::FSD::FSD
(
    const word& modelType,
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    reactionRateFlameArea_
    (
        reactionRateFlameArea::New
        (
            this->coeffs(),
            this->mesh(),
            *this
        )
    ),
    ft_
    (
        IOobject
        (
            this->thermo().phasePropertyName("ft"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    YFuelFuelStream_(dimensionedScalar(dimless, 1.0)),
    YO2OxiStream_(dimensionedScalar(dimless, 0.23)),
    Cv_(this->coeffs().template lookup<scalar>("Cv")),
    C_(5.0),
    ftMin_(0.0),
    ftMax_(1.0),
    ftDim_(300),
    ftVarMin_(this->coeffs().template lookup<scalar>("ftVarMin"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::FSD::~FSD()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::FSD::calculateSourceNorm()
{
    this->fresCorrect();

    const label fuelI = this->fuelIndex();

    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];

    const volScalarField& YO2 = this->thermo().composition().Y("O2");

    const dimensionedScalar s = this->s();

    ft_ =
        (s*YFuel - (YO2 - YO2OxiStream_))/(s*YFuelFuelStream_ + YO2OxiStream_);


    volVectorField nft(fvc::grad(ft_));

    volScalarField mgft(mag(nft));

    surfaceVectorField SfHat(this->mesh().Sf()/this->mesh().magSf());

    volScalarField cAux(scalar(1) - ft_);

    dimensionedScalar dMgft = 1.0e-3*
        (ft_*cAux*mgft)().weightedAverage(this->mesh().V())
       /((ft_*cAux)().weightedAverage(this->mesh().V()) + small)
      + dimensionedScalar(mgft.dimensions(), small);

    mgft += dMgft;

    nft /= mgft;

    const volVectorField& U = YO2.db().lookupObject<volVectorField>("U");

    const volScalarField sigma
    (
        (nft & nft)*fvc::div(U) - (nft & fvc::grad(U) & nft)
    );

    reactionRateFlameArea_->correct(sigma);

    const volScalarField& omegaFuel = reactionRateFlameArea_->omega();


    const scalar ftStoich =
        YO2OxiStream_.value()
       /(
            s.value()*YFuelFuelStream_.value() + YO2OxiStream_.value()
        );

    tmp<volScalarField> tPc
    (
        volScalarField::New
        (
            this->thermo().phasePropertyName("Pc"),
            U.mesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField& pc = tPc.ref();

    tmp<volScalarField> tomegaFuel
    (
        volScalarField::New
        (
            this->thermo().phasePropertyName("omegaFuelBar"),
            U.mesh(),
            dimensionedScalar(omegaFuel.dimensions(), 0)
        )
    );

    volScalarField& omegaFuelBar = tomegaFuel.ref();

    // Calculation of the mixture fraction variance (ftVar)
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
            momentumTransportModel::typeName
        );

    const volScalarField& delta = lesModel.delta();
    const volScalarField ftVar(Cv_*sqr(delta)*sqr(mgft));

    // Thickened flame (average flame thickness for counterflow configuration
    // is 1.5 mm)

    volScalarField  deltaF
    (
        lesModel.delta()/dimensionedScalar(dimLength, 1.5e-3)
    );

    // Linear correlation between delta and flame thickness
    volScalarField omegaF(max(deltaF*(4.0/3.0) + (2.0/3.0), scalar(1)));

    scalar deltaFt = 1.0/ftDim_;

    forAll(ft_, celli)
    {
        if (ft_[celli] > ftMin_ && ft_[celli] < ftMax_)
        {
            scalar ftCell = ft_[celli];

            if (ftVar[celli] > ftVarMin_) // sub-grid beta pdf of ft_
            {
                scalar ftVarc = ftVar[celli];
                scalar a =
                    max(ftCell*(ftCell*(1.0 - ftCell)/ftVarc - 1.0), 0.0);
                scalar b = max(a/ftCell - a, 0.0);

                for (int i=1; i<ftDim_; i++)
                {
                    scalar ft = i*deltaFt;
                    pc[celli] += pow(ft, a-1.0)*pow(1.0 - ft, b - 1.0)*deltaFt;
                }

                for (int i=1; i<ftDim_; i++)
                {
                    scalar ft = i*deltaFt;
                    omegaFuelBar[celli] +=
                        omegaFuel[celli]/omegaF[celli]
                       *exp
                        (
                           -sqr(ft - ftStoich)
                           /(2.0*sqr(0.01*omegaF[celli]))
                        )
                       *pow(ft, a - 1.0)
                       *pow(1.0 - ft, b - 1.0)
                       *deltaFt;
                }
                omegaFuelBar[celli] /= max(pc[celli], 1e-4);
            }
            else
            {
                omegaFuelBar[celli] =
                   omegaFuel[celli]/omegaF[celli]
                  *exp(-sqr(ftCell - ftStoich)/(2.0*sqr(0.01*omegaF[celli])));
            }
        }
        else
        {
            omegaFuelBar[celli] = 0.0;
        }
    }


    // Combustion progress variable, c

    List<label> productsIndex(2, label(-1));
    {
        label i = 0;
        forAll(this->specieProd(), specieI)
        {
            if (this->specieProd()[specieI] < 0)
            {
                productsIndex[i] = specieI;
                i++;
            }
        }
    }


    // Flamelet probability of the progress c based on IFC (reuse pc)
    scalar YprodTotal = 0;
    forAll(productsIndex, j)
    {
        YprodTotal += this->Yprod0()[productsIndex[j]];
    }

    forAll(ft_, celli)
    {
        if (ft_[celli] < ftStoich)
        {
            pc[celli] = ft_[celli]*(YprodTotal/ftStoich);
        }
        else
        {
            pc[celli] = (1.0 - ft_[celli])*(YprodTotal/(1.0 - ftStoich));
        }
    }

    tmp<volScalarField> tproducts
    (
        volScalarField::New
        (
            this->thermo().phasePropertyName("products"),
            U.mesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField& products = tproducts.ref();

    forAll(productsIndex, j)
    {
        label specieI = productsIndex[j];
        const volScalarField& Yp = this->thermo().composition().Y()[specieI];
        products += Yp;
    }

    volScalarField c
    (
        max(scalar(1) - products/max(pc, scalar(1e-5)), scalar(0))
    );

    pc = min(C_*c, scalar(1));

    const volScalarField fres(this->fres(fuelI));

    this->wFuel_ == mgft*pc*omegaFuelBar;
}


void Foam::combustionModels::FSD::correct()
{
    this->wFuel_ ==
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0);

    calculateSourceNorm();
}


bool Foam::combustionModels::FSD::read()
{
    if (singleStepCombustion::read())
    {
        this->coeffs().lookup("Cv") >> Cv_ ;
        this->coeffs().lookup("ftVarMin") >> ftVarMin_;
        reactionRateFlameArea_->read(this->coeffs());
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

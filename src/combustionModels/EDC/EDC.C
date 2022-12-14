/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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

#include "EDC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum<Foam::combustionModels::EDCversions, 4>::names[] =
    {"v1981", "v1996", "v2005", "v2016"};

const Foam::NamedEnum<Foam::combustionModels::EDCversions, 4>
    Foam::combustionModels::EDCversionNames;

const Foam::combustionModels::EDCversions
    Foam::combustionModels::EDCdefaultVersion =
    Foam::combustionModels::EDCversions::v2005;

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(EDC, 0);
    addToRunTimeSelectionTable(combustionModel, EDC, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::EDC::EDC
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb, combustionProperties),
    version_
    (
        EDCversionNames
        [
            this->coeffs().lookupOrDefault
            (
                "version",
                word(EDCversionNames[EDCdefaultVersion])
            )
        ]
    ),
    C1_(this->coeffs().lookupOrDefault("C1", 0.05774)),
    C2_(this->coeffs().lookupOrDefault("C2", 0.5)),
    Cgamma_(this->coeffs().lookupOrDefault("Cgamma", 2.1377)),
    Ctau_(this->coeffs().lookupOrDefault("Ctau", 0.4083)),
    exp1_(this->coeffs().lookupOrDefault("exp1", EDCexp1[int(version_)])),
    exp2_(this->coeffs().lookupOrDefault("exp2", EDCexp2[int(version_)])),
    kappa_
    (
        IOobject
        (
            this->thermo().phasePropertyName(typedName("kappa")),
            this->mesh().time().name(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    outerCorrect_
    (
        this->coeffs().lookupOrDefault("outerCorrect", true)
    ),
    timeIndex_(-1),
    chemistryPtr_(basicChemistryModel::New(thermo))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::EDC::~EDC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::EDC::correct()
{
    if (!outerCorrect_ && timeIndex_ == this->mesh().time().timeIndex())
    {
        return;
    }

    tmp<volScalarField> tepsilon(this->turbulence().epsilon());
    const volScalarField& epsilon = tepsilon();

    tmp<volScalarField> tnu(this->turbulence().nu());
    const volScalarField& nu = tnu();

    tmp<volScalarField> tk(this->turbulence().k());
    const volScalarField& k = tk();

    scalarField tauStar(epsilon.size(), 0);

    if (version_ == EDCversions::v2016)
    {
        tmp<volScalarField> ttc(chemistryPtr_->tc());
        const volScalarField& tc = ttc();

        forAll(tauStar, i)
        {
            const scalar Da =
                max(min(sqrt(nu[i]/(epsilon[i] + small))/tc[i], 10), 1e-10);

            const scalar ReT = sqr(k[i])/(nu[i]*epsilon[i] + small);
            const scalar CtauI = min(C1_/(Da*sqrt(ReT + 1)), Ctau_);

            const scalar CgammaI =
                max(min(C2_*sqrt(Da*(ReT + 1)), 5), Cgamma_);

            const scalar gammaL =
                CgammaI*pow025(nu[i]*epsilon[i]/(sqr(k[i]) + small));

            tauStar[i] = CtauI*sqrt(nu[i]/(epsilon[i] + small));

            if (gammaL >= 1)
            {
                kappa_[i] = 1;
            }
            else
            {
                kappa_[i] =
                    max
                    (
                        min
                        (
                            pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                            1
                        ),
                        0
                    );
            }
        }
    }
    else
    {
        forAll(tauStar, i)
        {
            const scalar gammaL =
                Cgamma_*pow025(nu[i]*epsilon[i]/(sqr(k[i]) + small));

            tauStar[i] = Ctau_*sqrt(nu[i]/(epsilon[i] + small));
            if (gammaL >= 1)
            {
                kappa_[i] = 1;
            }
            else
            {
                kappa_[i] =
                    max
                    (
                        min
                        (
                            pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                            1
                        ),
                        0
                    );
            }
        }
    }

    chemistryPtr_->solve(tauStar);

    timeIndex_ = this->mesh().time().timeIndex();
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::EDC::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    const label specieI = this->thermo().composition().species()[Y.member()];
    Su += chemistryPtr_->RR()[specieI];

    return kappa_*tSu;
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::EDC::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typedName("Qdot")),
        kappa_*chemistryPtr_->Qdot()
    );
}


bool Foam::combustionModels::EDC::read()
{
    if (combustionModel::read())
    {
        version_ =
        (
            EDCversionNames
            [
                this->coeffs().lookupOrDefault
                (
                    "version",
                    word(EDCversionNames[EDCdefaultVersion])
                )
            ]
        );
        C1_ = this->coeffs().lookupOrDefault("C1", 0.05774);
        C2_ = this->coeffs().lookupOrDefault("C2", 0.5);
        Cgamma_ = this->coeffs().lookupOrDefault("Cgamma", 2.1377);
        Ctau_ = this->coeffs().lookupOrDefault("Ctau", 0.4083);
        exp1_ = this->coeffs().lookupOrDefault("exp1", EDCexp1[int(version_)]);
        exp2_ = this->coeffs().lookupOrDefault("exp2", EDCexp2[int(version_)]);
        outerCorrect_ = this->coeffs().lookupOrDefault("outerCorrect", true);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

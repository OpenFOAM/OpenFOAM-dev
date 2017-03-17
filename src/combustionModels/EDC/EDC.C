/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::EDC<Type>::EDC
(
    const word& modelType,
    const fvMesh& mesh,
    const word& combustionProperties,
    const word& phaseName
)
:
    laminar<Type>(modelType, mesh, combustionProperties, phaseName),
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
            IOobject::groupName(typeName + ":kappa", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::EDC<Type>::~EDC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::combustionModels::EDC<Type>::correct()
{
    if (this->active())
    {
        tmp<volScalarField> tepsilon(this->turbulence().epsilon());
        const volScalarField& epsilon = tepsilon();

        tmp<volScalarField> tmu(this->turbulence().mu());
        const volScalarField& mu = tmu();

        tmp<volScalarField> tk(this->turbulence().k());
        const volScalarField& k = tk();

        tmp<volScalarField> trho(this->rho());
        const volScalarField& rho = trho();

        scalarField tauStar(epsilon.size(), 0);

        if (version_ == EDCversions::v2016)
        {
            tmp<volScalarField> ttc(this->chemistryPtr_->tc());
            const volScalarField& tc = ttc();

            forAll(tauStar, i)
            {
                const scalar nu = mu[i]/(rho[i] + SMALL);

                const scalar Da =
                    max(min(sqrt(nu/(epsilon[i] + SMALL))/tc[i], 10), 1e-10);

                const scalar ReT = sqr(k[i])/(nu*epsilon[i] + SMALL);
                const scalar CtauI = min(C1_/(Da*sqrt(ReT + 1)), 2.1377);

                const scalar CgammaI =
                    max(min(C2_*sqrt(Da*(ReT + 1)), 5), 0.4082);

                const scalar gammaL =
                    CgammaI*pow025(nu*epsilon[i]/(sqr(k[i]) + SMALL));

                tauStar[i] = CtauI*sqrt(nu/(epsilon[i] + SMALL));

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
                const scalar nu = mu[i]/(rho[i] + SMALL);
                const scalar gammaL =
                    Cgamma_*pow025(nu*epsilon[i]/(sqr(k[i]) + SMALL));

                tauStar[i] = Ctau_*sqrt(nu/(epsilon[i] + SMALL));
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

        this->chemistryPtr_->solve(tauStar);
    }
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::EDC<Type>::R(volScalarField& Y) const
{
    return kappa_*laminar<Type>::R(Y);
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::EDC<Type>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName(typeName + ":Qdot", this->phaseName_),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->active())
    {
        tQdot.ref() = kappa_*this->chemistryPtr_->Qdot();
    }

    return tQdot;
}


template<class Type>
bool Foam::combustionModels::EDC<Type>::read()
{
    if (Type::read())
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

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

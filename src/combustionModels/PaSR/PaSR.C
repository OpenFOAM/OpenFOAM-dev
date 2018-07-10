/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "PaSR.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::PaSR<ReactionThermo>::PaSR
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    Cmix_(readScalar(this->coeffs().lookup("Cmix"))),
    kappa_
    (
        IOobject
        (
            thermo.phasePropertyName(typeName + ":kappa"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("kappa", dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::PaSR<ReactionThermo>::~PaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::PaSR<ReactionThermo>::correct()
{
    if (this->active())
    {
        laminar<ReactionThermo>::correct();

        tmp<volScalarField> tepsilon(this->turbulence().epsilon());
        const scalarField& epsilon = tepsilon();

        tmp<volScalarField> tmuEff(this->turbulence().muEff());
        const scalarField& muEff = tmuEff();

        tmp<volScalarField> ttc(this->tc());
        const scalarField& tc = ttc();

        tmp<volScalarField> trho(this->rho());
        const scalarField& rho = trho();

        forAll(epsilon, i)
        {
            const scalar tk =
                Cmix_*sqrt(max(muEff[i]/rho[i]/(epsilon[i] + small), 0));

            if (tk > small)
            {
                kappa_[i] = tc[i]/(tc[i] + tk);
            }
            else
            {
                kappa_[i] = 1.0;
            }
        }
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR<ReactionThermo>::R(volScalarField& Y) const
{
    return kappa_*laminar<ReactionThermo>::R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            this->thermo().phasePropertyName(typeName + ":Qdot"),
            kappa_*laminar<ReactionThermo>::Qdot()
        )
    );
}


template<class ReactionThermo>
bool Foam::combustionModels::PaSR<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        this->coeffs().lookup("Cmix") >> Cmix_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

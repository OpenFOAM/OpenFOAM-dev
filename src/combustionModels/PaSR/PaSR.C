/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

template<class Type>
Foam::combustionModels::PaSR<Type>::PaSR
(
    const word& modelType,
    const fvMesh& mesh,
    const word& combustionProperties,
    const word& phaseName
)
:
    laminar<Type>(modelType, mesh, combustionProperties, phaseName),
    Cmix_(readScalar(this->coeffs().lookup("Cmix"))),
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
Foam::combustionModels::PaSR<Type>::~PaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::combustionModels::PaSR<Type>::correct()
{
    if (this->active())
    {
        laminar<Type>::correct();

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
                Cmix_*sqrt(max(muEff[i]/rho[i]/(epsilon[i] + SMALL), 0));

            if (tk > SMALL)
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


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR<Type>::R(volScalarField& Y) const
{
    return kappa_*laminar<Type>::R(Y);
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<Type>::Qdot() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName(typeName + ":Qdot", this->phaseName_),
            kappa_*laminar<Type>::Qdot()
        )
    );
}


template<class Type>
bool Foam::combustionModels::PaSR<Type>::read()
{
    if (laminar<Type>::read())
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

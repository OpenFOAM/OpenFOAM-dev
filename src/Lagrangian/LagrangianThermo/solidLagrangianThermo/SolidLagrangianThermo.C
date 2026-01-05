/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "SolidLagrangianThermo.H"
#include "uniformGeometricFields.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::SolidLagrangianThermo<BaseThermo>::~SolidLagrangianThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
void Foam::SolidLagrangianThermo<BaseThermo>::correct
(
    const LagrangianSubMesh& subMesh
)
{
    if (BaseThermo::debug) InfoInFunction << endl;

    const SubField<scalar> e = subMesh.sub(this->e_.primitiveField());

    const UniformField<scalar>& p = this->p_.primitiveField();
    SubField<scalar> T = subMesh.sub(this->T_.primitiveFieldRef());
    SubField<scalar> rho = subMesh.sub(this->rho_.primitiveFieldRef());
    SubField<scalar> Cv = subMesh.sub(this->Cv_.primitiveFieldRef());
    SubField<scalar> kappa = subMesh.sub(this->kappa_.primitiveFieldRef());

    auto Yslicer = this->Yslicer();

    forAll(T, subi)
    {
        const label i = subMesh.start() + subi;

        auto composition = this->elementComposition(Yslicer, i);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture =
            this->transportMixture(composition, thermoMixture);

        T[subi] = thermoMixture.Tes(e[subi], p[subi], T[subi]);

        rho[subi] = thermoMixture.rho(p[subi], T[subi]);
        Cv[subi] = thermoMixture.Cv(p[subi], T[subi]);

        kappa[subi] = transportMixture.kappa(p[subi], T[subi]);
    }

    if (BaseThermo::debug) Info<< "    Finished" << endl;
}


// ************************************************************************* //

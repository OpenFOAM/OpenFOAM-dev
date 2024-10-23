/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

#include "phasePressureModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::phasePressureModel::phasePressureModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<phaseCompressible::momentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    phase_(refCast<const phaseModel>(viscosity)),

    preAlphaExp_(coeffDict().lookup<scalar>("preAlphaExp")),
    expMax_(coeffDict().lookup<scalar>("expMax")),
    g0_
    (
        "g0",
        dimensionSet(1, -1, -2, 0, 0),
        coeffDict().lookup("g0")
    )
{
    nut_ == dimensionedScalar(nut_.dimensions(), 0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::phasePressureModel::~phasePressureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::phasePressureModel::read()
{
    if
    (
        eddyViscosity<RASModel<phaseCompressible::momentumTransportModel>>::
        read()
    )
    {
        coeffDict().lookup("preAlphaExp") >> preAlphaExp_;
        coeffDict().lookup("expMax") >> expMax_;
        g0_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::omega() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::R() const
{
    return tmp<volSymmTensorField>
    (
        volSymmTensorField::New
        (
            IOobject::groupName("R", U_.group()),
            mesh_,
            dimensioned<symmTensor>
            (
                "R",
                dimensionSet(0, 2, -2, 0, 0),
                Zero
            )
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::pPrime() const
{
    tmp<volScalarField> tpPrime
    (
        volScalarField::New
        (
            IOobject::groupName("pPrime", U_.group()),
            g0_
           *min
            (
                exp(preAlphaExp_*(alpha_ - phase_.alphaMax())),
                expMax_
            )
        )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


// Foam::tmp<Foam::surfaceScalarField>
// Foam::RASModels::phasePressureModel::pPrimef() const
// {
//     tmp<surfaceScalarField> tpPrime
//     (
//         surfaceScalarField::New
//         (
//             IOobject::groupName("pPrimef", U_.group()),
//             g0_
//            *min
//             (
//                 exp(preAlphaExp_
//                 *(fvc::interpolate(alpha_, "hmm") - phase_.alphaMax())),
//                 expMax_
//             )
//         )
//     );

//    surfaceScalarField::Boundary& bpPrime =
//        tpPrime.ref().boundaryFieldRef();

//     forAll(bpPrime, patchi)
//     {
//         if (!bpPrime[patchi].coupled())
//         {
//             bpPrime[patchi] == 0;
//         }
//     }

//     return tpPrime;
// }


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::phasePressureModel::pPrimef() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("pPrimef", U_.group()),
        fvc::interpolate(pPrime())
    );
}


Foam::tmp<Foam::surfaceVectorField>
Foam::RASModels::phasePressureModel::devTau() const
{
    return tmp<surfaceVectorField>
    (
        surfaceVectorField::New
        (
            IOobject::groupName("devTau", U_.group()),
            mesh_,
            dimensioned<vector>
            (
                "devTau",
                rho_.dimensions()*dimensionSet(0, 2, -2, 0, 0),
                Zero
            )
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::phasePressureModel::divDevTau
(
    volVectorField& U
) const
{
    return tmp<fvVectorMatrix>
    (
        new fvVectorMatrix
        (
            U,
            rho_.dimensions()*dimensionSet(0, 4, -2, 0, 0)
        )
    );
}


void Foam::RASModels::phasePressureModel::correct()
{}


// ************************************************************************* //

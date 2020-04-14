/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "eddyViscosity.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::eddyViscosity<BasicMomentumTransportModel>::eddyViscosity
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
:
    linearViscousStress<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    nut_
    (
        IOobject
        (
            IOobject::groupName("nut", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool Foam::eddyViscosity<BasicMomentumTransportModel>::read()
{
    return BasicMomentumTransportModel::read();
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::eddyViscosity<BasicMomentumTransportModel>::sigma() const
{
    tmp<volScalarField> tk(k());

    // Get list of patchField type names from k
    wordList patchFieldTypes(tk().boundaryField().types());

    // For k patchField types which do not have an equivalent for symmTensor
    // set to calculated
    forAll(patchFieldTypes, i)
    {
        if
        (
           !fvPatchField<symmTensor>::patchConstructorTablePtr_
                ->found(patchFieldTypes[i])
        )
        {
            patchFieldTypes[i] = calculatedFvPatchField<symmTensor>::typeName;
        }
    }

    return volSymmTensorField::New
    (
        IOobject::groupName("R", this->alphaRhoPhi_.group()),
        ((2.0/3.0)*I)*tk() - (nut_)*dev(twoSymm(fvc::grad(this->U_))),
        patchFieldTypes
    );
}


template<class BasicMomentumTransportModel>
void Foam::eddyViscosity<BasicMomentumTransportModel>::validate()
{
    correctNut();
}


template<class BasicMomentumTransportModel>
void Foam::eddyViscosity<BasicMomentumTransportModel>::correct()
{
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //

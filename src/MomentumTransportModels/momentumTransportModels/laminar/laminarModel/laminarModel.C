/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "laminarModel.H"
#include "Stokes.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void Foam::laminarModel<BasicMomentumTransportModel>::printCoeffs
(
    const word& type
)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::laminarModel<BasicMomentumTransportModel>::laminarModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    BasicMomentumTransportModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    laminarDict_(this->subOrEmptyDict("laminar")),
    printCoeffs_(laminarDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(laminarDict_.optionalSubDict(type + "Coeffs"))
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::autoPtr<Foam::laminarModel<BasicMomentumTransportModel>>
Foam::laminarModel<BasicMomentumTransportModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
{
    const IOdictionary modelDict
    (
        momentumTransportModel::readModelDict
        (
            U.db(),
            alphaRhoPhi.group()
        )
    );

    if (modelDict.found("laminar"))
    {
        const word modelType =
            modelDict.subDict("laminar").lookupBackwardsCompatible<word>
            (
                {"model", "laminarModel"}
            );

        Info<< "Selecting laminar stress model " << modelType << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown laminarModel type "
                << modelType << nl << nl
                << "Valid laminarModel types:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<laminarModel>
        (
            cstrIter()
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                viscosity
            )
        );
    }
    else
    {
        Info<< "Selecting laminar stress model "
            << laminarModels::Stokes<BasicMomentumTransportModel>::typeName
            << endl;

        return autoPtr<laminarModel>
        (
            new laminarModels::Stokes<BasicMomentumTransportModel>
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                viscosity
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool Foam::laminarModel<BasicMomentumTransportModel>::read()
{
    if (BasicMomentumTransportModel::read())
    {
        laminarDict_ <<= this->subDict("laminar");

        coeffDict_ <<= laminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMomentumTransportModel>::nut() const
{
    return volScalarField::New
    (
        this->groupName("nut"),
        this->mesh_,
        dimensionedScalar(dimViscosity, 0)
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::scalarField>
Foam::laminarModel<BasicMomentumTransportModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), 0.0)
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMomentumTransportModel>::k() const
{
    return volScalarField::New
    (
        this->groupName("k"),
        this->mesh_,
        dimensionedScalar(sqr(dimVelocity), 0)
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMomentumTransportModel>::epsilon() const
{
    return volScalarField::New
    (
        this->groupName("epsilon"),
        this->mesh_,
        dimensionedScalar(sqr(dimVelocity)/dimTime, 0)
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMomentumTransportModel>::omega() const
{
    return volScalarField::New
    (
        this->groupName("omega"),
        this->mesh_,
        dimensionedScalar(dimless/dimTime, 0)
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::laminarModel<BasicMomentumTransportModel>::sigma() const
{
    return volSymmTensorField::New
    (
        this->groupName("sigma"),
        this->mesh_,
        dimensionedSymmTensor(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicMomentumTransportModel>
void Foam::laminarModel<BasicMomentumTransportModel>::predict()
{
    BasicMomentumTransportModel::predict();
}


template<class BasicMomentumTransportModel>
void Foam::laminarModel<BasicMomentumTransportModel>::correct()
{
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //

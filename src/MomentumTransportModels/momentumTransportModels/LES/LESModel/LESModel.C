/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "LESModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void Foam::LESModel<BasicMomentumTransportModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::LESModel<BasicMomentumTransportModel>::LESModel
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
    BasicMomentumTransportModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    LESDict_(this->subOrEmptyDict("LES")),
    turbulence_(LESDict_.lookup("turbulence")),
    printCoeffs_(LESDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(LESDict_.optionalSubDict(type + "Coeffs")),

    kMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kMin",
            LESDict_,
            sqr(dimVelocity),
            small
        )
    ),

    epsilonMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "epsilonMin",
            LESDict_,
            kMin_.dimensions()/dimTime,
            small
        )
    ),

    omegaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaMin",
            LESDict_,
            dimless/dimTime,
            small
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            LESDict_
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::autoPtr<Foam::LESModel<BasicMomentumTransportModel>>
Foam::LESModel<BasicMomentumTransportModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
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

    const word modelType =
        modelDict.subDict("LES").lookupBackwardsCompatible<word>
        (
            {"model", "LESModel"}
        );

    Info<< "Selecting LES turbulence model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown LESModel type "
            << modelType << nl << nl
            << "Valid LESModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<LESModel>
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, transport)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool Foam::LESModel<BasicMomentumTransportModel>::read()
{
    if (BasicMomentumTransportModel::read())
    {
        LESDict_ <<= this->subDict("LES");
        LESDict_.lookup("turbulence") >> turbulence_;

        coeffDict_ <<= LESDict_.optionalSubDict(type() + "Coeffs");

        delta_().read(LESDict_);

        kMin_.readIfPresent(LESDict_);

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void Foam::LESModel<BasicMomentumTransportModel>::correct()
{
    delta_().correct();
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //

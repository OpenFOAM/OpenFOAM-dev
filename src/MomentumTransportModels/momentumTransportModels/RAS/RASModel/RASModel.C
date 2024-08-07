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

#include "RASModel.H"
#include "NewtonianViscosityModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::RASModel<BasicMomentumTransportModel>::RASModel
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

    turbulence_(RASDict().lookup("turbulence")),
    kMin_("kMin", sqr(dimVelocity), RASDict(), small),
    nutMaxCoeff_("nutMaxCoeff", dimless, RASDict(), 1e5),

    viscosityModel_
    (
        coeffDict().found("viscosityModel")
      ? laminarModels::generalisedNewtonianViscosityModel::New
        (
            coeffDict(),
            viscosity,
            U
        )
      : autoPtr<laminarModels::generalisedNewtonianViscosityModel>
        (
            new laminarModels::generalisedNewtonianViscosityModels::Newtonian
            (
                coeffDict(),
                viscosity,
                U
            )
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::autoPtr<Foam::RASModel<BasicMomentumTransportModel>>
Foam::RASModel<BasicMomentumTransportModel>::New
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

    const word modelType =
        modelDict.subDict("RAS").lookupBackwardsCompatible<word>
        (
            {"model", "RASModel"}
        );

    Info<< indent
        << "Selecting RAS turbulence model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown RASModel type "
            << modelType << nl << nl
            << "Valid RASModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    Info<< incrIndent;

    autoPtr<RASModel> modelPtr
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, viscosity)
    );

    Info<< decrIndent;

    return modelPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
const Foam::dictionary&
Foam::RASModel<BasicMomentumTransportModel>::RASDict() const
{
    return this->subDict("RAS");
}


template<class BasicMomentumTransportModel>
const Foam::dictionary&
Foam::RASModel<BasicMomentumTransportModel>::coeffDict() const
{
    return this->RASDict().optionalSubDict(type() + "Coeffs");
}


template<class BasicMomentumTransportModel>
bool Foam::RASModel<BasicMomentumTransportModel>::read()
{
    if (BasicMomentumTransportModel::read())
    {
        RASDict().lookup("turbulence") >> turbulence_;
        kMin_.readIfPresent(RASDict());
        nutMaxCoeff_.readIfPresent(RASDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void Foam::RASModel<BasicMomentumTransportModel>::correct()
{
    viscosityModel_->correct();
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //

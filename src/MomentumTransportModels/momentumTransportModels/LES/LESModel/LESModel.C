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

#include "LESModel.H"
#include "NewtonianViscosityModel.H"

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

    turbulence_(LESDict().lookup("turbulence")),
    kMin_("kMin", sqr(dimVelocity), LESDict(), small),
    nutMaxCoeff_("nutMaxCoeff", dimless, LESDict(), 1e5),

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
    ),

    delta_
    (
        LESdelta::New
        (
            this->groupName("delta"),
            *this,
            LESDict()
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
        modelDict.subDict("LES").lookupBackwardsCompatible<word>
        (
            {"model", "LESModel"}
        );

    Info<< indent
        << "Selecting LES turbulence model " << modelType << endl;

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

    Info<< incrIndent;

    autoPtr<LESModel> modelPtr
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, viscosity)
    );

    Info<< decrIndent;

    return modelPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
const Foam::dictionary&
Foam::LESModel<BasicMomentumTransportModel>::LESDict() const
{
    return this->subDict("LES");
}


template<class BasicMomentumTransportModel>
const Foam::dictionary&
Foam::LESModel<BasicMomentumTransportModel>::coeffDict() const
{
    return this->LESDict().optionalSubDict(type() + "Coeffs");
}


template<class BasicMomentumTransportModel>
bool Foam::LESModel<BasicMomentumTransportModel>::read()
{
    if (BasicMomentumTransportModel::read())
    {
        LESDict().lookup("turbulence") >> turbulence_;
        delta_().read(LESDict());

        kMin_.readIfPresent(LESDict());
        nutMaxCoeff_.readIfPresent(LESDict());

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
    viscosityModel_->correct();
    delta_().correct();
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //

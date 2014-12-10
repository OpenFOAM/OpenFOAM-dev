/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< type << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::LESModel<BasicTurbulenceModel>::LESModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    LESDict_(this->subOrEmptyDict("LES")),
    turbulence_(LESDict_.lookup("turbulence")),
    printCoeffs_(LESDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(LESDict_.subOrEmptyDict(type + "Coeffs")),

    kMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kMin",
            LESDict_,
            SMALL,
            sqr(dimVelocity)
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", U.group()),
            U.mesh(),
            LESDict_
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::LESModel<BasicTurbulenceModel> >
Foam::LESModel<BasicTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                IOobject::groupName(propertiesName, U.group()),
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("LES").lookup("LESModel")
    );

    Info<< "Selecting LES turbulence model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LESModel::New"
            "("
                "const volScalarField&, "
                "const volVectorField&, "
                "const surfaceScalarField&, "
                "transportModel&, "
                "const word&"
            ")"
        )   << "Unknown LESModel type "
            << modelType << nl << nl
            << "Valid LESModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<LESModel>
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


template<class BasicTurbulenceModel>
bool Foam::LESModel<BasicTurbulenceModel>::read()
{
    if (turbulenceModel::read())
    {
        LESDict_ <<= this->subDict("LES");
        LESDict_.lookup("turbulence") >> turbulence_;

        if (const dictionary* dictPtr = LESDict_.subDictPtr(type() + "Coeffs"))
        {
            coeffDict_ <<= *dictPtr;
        }

        kMin_.readIfPresent(LESDict_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

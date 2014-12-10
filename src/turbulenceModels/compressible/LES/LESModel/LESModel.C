/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LESModel, 0);
defineRunTimeSelectionTable(LESModel, dictionary);
addToRunTimeSelectionTable(turbulenceModel, LESModel, turbulenceModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void LESModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

LESModel::LESModel
(
    const word& type,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermoPhysicalModel,
    const word& turbulenceModelName
)
:
    turbulenceModel(rho, U, phi, thermoPhysicalModel, turbulenceModelName),

    IOdictionary
    (
        IOobject
        (
            "LESProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subOrEmptyDict(type + "Coeffs")),

    kMin_("kMin", sqr(dimVelocity), SMALL),

    delta_(LESdelta::New("delta", U.mesh(), *this))
{
    kMin_.readIfPresent(*this);

    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<LESModel> LESModel::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermoPhysicalModel,
    const word& turbulenceModelName
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
                "LESProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("LESModel")
    );

    Info<< "Selecting LES turbulence model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
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
                "const fluidThermo&, "
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
        cstrIter()(rho, U, phi, thermoPhysicalModel, turbulenceModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LESModel::correct(const tmp<volTensorField>&)
{
    turbulenceModel::correct();
    delta_().correct();
}


void LESModel::correct()
{
    correct(fvc::grad(U_));
}


bool LESModel::read()
{
    // Bit of trickery : we are both IOdictionary ('RASProperties') and
    // an regIOobject (from the turbulenceModel). Problem is to distinguish
    // between the two - we only want to reread the IOdictionary.

    bool ok = IOdictionary::readData
    (
        IOdictionary::readStream
        (
            IOdictionary::type()
        )
    );
    IOdictionary::close();

    if (ok)
    {
        if (const dictionary* dictPtr = subDictPtr(type() + "Coeffs"))
        {
            coeffDict_ <<= *dictPtr;
        }

        kMin_.readIfPresent(*this);

        delta_().read(*this);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //

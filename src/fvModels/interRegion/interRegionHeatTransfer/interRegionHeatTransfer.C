/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "interRegionHeatTransfer.H"
#include "basicThermo.H"
#include "fvmSup.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcVolumeIntegrate.H"
#include "fvModels.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        interRegionHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::interRegionHeatTransfer::readCoeffs()
{
    semiImplicit_ = coeffs().lookup<bool>("semiImplicit");

    TName_ = coeffs().lookupOrDefault<word>("T", "T");
    TNbrName_ = coeffs().lookupOrDefault<word>("TNbr", "T");

    if (master())
    {
        heatTransferModel_ = heatTransferModel::New(coeffs(), *this);
    }
}


const Foam::fv::heatTransferModel&
Foam::fv::interRegionHeatTransfer::nbrHeatTransferModel() const
{
    return
        refCast<const interRegionHeatTransfer>(nbrModel()).heatTransferModel_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransfer::interRegionHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionModel(name, modelType, dict, mesh),
    semiImplicit_(false),
    TName_(word::null),
    TNbrName_(word::null),
    heatTransferModel_(nullptr)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransfer::~interRegionHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::interRegionHeatTransfer::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::interRegionHeatTransfer::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const volScalarField& he = eqn.psi();

    const volScalarField& T =
        mesh().lookupObject<volScalarField>(TName_);
    const volScalarField& Tnbr =
        nbrMesh().lookupObject<volScalarField>(TNbrName_);

    tmp<volScalarField> tTnbrMapped =
        volScalarField::New(TName_ + "nbrMapped", T);
    interpolate(Tnbr, tTnbrMapped->primitiveFieldRef());
    volScalarField& TnbrMapped = tTnbrMapped.ref();

    // Get the heat transfer coefficient field
    tmp<volScalarField> tHtcAoV;
    if (master())
    {
        tmp<volScalarField> mask =
            volScalarField::New
            (
                "mask",
                mesh(),
                dimensionedScalar(dimless, 0)
            );
        tmp<volScalarField> oneNbr =
            volScalarField::New
            (
                "one",
                nbrMesh(),
                dimensionedScalar(dimless, 1)
            );
        interpolate(oneNbr(), mask.ref().primitiveFieldRef());
        tHtcAoV =
            mask
           *heatTransferModel_->htc()
           *heatTransferModel_->AoV();
    }
    else
    {
        tmp<volScalarField> tHtcNbr =
            nbrHeatTransferModel().htc()
           *nbrHeatTransferModel().AoV();
        tHtcAoV =
            volScalarField::New
            (
                tHtcNbr().name(),
                mesh(),
                dimensionedScalar(tHtcNbr().dimensions(), 0)
            );
        interpolate(tHtcNbr(), tHtcAoV.ref().primitiveFieldRef());
    }
    const volScalarField& htcAoV = tHtcAoV();

    if (semiImplicit_)
    {
        if (he.dimensions() == dimEnergy/dimMass)
        {
            const basicThermo& thermo =
               mesh().lookupObject<basicThermo>(physicalProperties::typeName);

            const volScalarField htcAoVByCpv(htcAoV/thermo.Cpv());

            eqn +=
                htcAoV*(TnbrMapped - T)
              + htcAoVByCpv*he - fvm::Sp(htcAoVByCpv, he);
        }
        else if (he.dimensions() == dimTemperature)
        {
            eqn += htcAoV*TnbrMapped - fvm::Sp(htcAoV, he);
        }
    }
    else
    {
        eqn += htcAoV*(TnbrMapped - T);
    }
}


void Foam::fv::interRegionHeatTransfer::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    addSup(eqn, fieldName);
}


void Foam::fv::interRegionHeatTransfer::correct()
{
    if (master())
    {
        heatTransferModel_->correct();
    }
}


bool Foam::fv::interRegionHeatTransfer::read(const dictionary& dict)
{
    if (interRegionModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

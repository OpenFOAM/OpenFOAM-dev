/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

void Foam::fv::interRegionHeatTransfer::readCoeffs(const dictionary& dict)
{
    semiImplicit_ = dict.lookup<bool>("semiImplicit");

    TName_ = dict.lookupOrDefault<word>("T", "T");
    TNbrName_ = dict.lookupOrDefault<word>("TNbr", "T");

    if (master())
    {
        heatTransferAv_.reset(new heatTransferAv(dict, mesh()));

        heatTransferCoefficientModel_ =
            heatTransferCoefficientModel::New(dict, *this);
    }
}


const Foam::fv::interRegionHeatTransfer&
Foam::fv::interRegionHeatTransfer::nbrHeatTransfer() const
{
    return refCast<const interRegionHeatTransfer>(nbrModel());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransfer::interRegionHeatTransfer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    interRegionModel(name, modelType, mesh, dict),
    semiImplicit_(false),
    TName_(word::null),
    TNbrName_(word::null),
    heatTransferAv_(nullptr),
    heatTransferCoefficientModel_(nullptr)
{
    readCoeffs(coeffs(dict));
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
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField& T =
        mesh().lookupObject<volScalarField>(TName_);

    tmp<volScalarField> tTnbr = volScalarField::New(TNbrName_, T);
    interpolate
    (
        nbrMesh().lookupObject<volScalarField>(TNbrName_),
        tTnbr->primitiveFieldRef()
    );
    const volScalarField& Tnbr = tTnbr();

    // Get the heat transfer coefficient field
    tmp<volScalarField> tHtcAv;
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
        tHtcAv =
            mask
           *heatTransferCoefficientModel_->htc()
           *heatTransferAv_->Av();
    }
    else
    {
        tmp<volScalarField> tHtcNbr =
            nbrHeatTransfer().heatTransferCoefficientModel_->htc()
           *nbrHeatTransfer().heatTransferAv_->Av();
        tHtcAv =
            volScalarField::New
            (
                tHtcNbr().name(),
                mesh(),
                dimensionedScalar(tHtcNbr().dimensions(), 0)
            );
        interpolate(tHtcNbr(), tHtcAv.ref().primitiveFieldRef());
    }
    const volScalarField& htcAv = tHtcAv();

    if (semiImplicit_)
    {
        if (he.dimensions() == dimEnergy/dimMass)
        {
            const basicThermo& thermo =
               mesh().lookupObject<basicThermo>(physicalProperties::typeName);

            const volScalarField htcAvByCpv(htcAv/thermo.Cpv());

            eqn +=
                htcAv*(Tnbr - T)
              + htcAvByCpv*he - fvm::Sp(htcAvByCpv, he);
        }
        else if (he.dimensions() == dimTemperature)
        {
            eqn += htcAv*Tnbr - fvm::Sp(htcAv, he);
        }
    }
    else
    {
        eqn += htcAv*(Tnbr - T);
    }
}


void Foam::fv::interRegionHeatTransfer::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    addSup(he, eqn);
}


void Foam::fv::interRegionHeatTransfer::correct()
{
    if (master())
    {
        heatTransferCoefficientModel_->correct();
    }
}


bool Foam::fv::interRegionHeatTransfer::movePoints()
{
    NotImplemented;
    return true;
}


void Foam::fv::interRegionHeatTransfer::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fv::interRegionHeatTransfer::mapMesh(const polyMeshMap&)
{
    NotImplemented;
}


void Foam::fv::interRegionHeatTransfer::distribute
(
    const polyDistributionMap&
)
{
    NotImplemented;
}


bool Foam::fv::interRegionHeatTransfer::read(const dictionary& dict)
{
    if (interRegionModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

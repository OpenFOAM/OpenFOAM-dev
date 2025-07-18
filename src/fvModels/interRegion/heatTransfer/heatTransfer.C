/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "heatTransfer.H"
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
    defineTypeNameAndDebug(heatTransfer, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        heatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransfer::readCoeffs(const dictionary& dict)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    semiImplicit_ = dict.lookup<bool>("semiImplicit");

    TName_ = dict.lookupOrDefault<word>("T", "T");

    Ta_ = dimensionedScalar("Ta", dimTemperature, dict);

    heatTransferAv_.reset(new heatTransferAv(dict, mesh()));

    heatTransferCoefficientModel_ =
        heatTransferCoefficientModel::New(dict, mesh());
}


template<class AlphaFieldType>
void Foam::fv::heatTransfer::add
(
    const AlphaFieldType& alpha,
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField& he = eqn.psi();

    const volScalarField& T =
        mesh().lookupObject<volScalarField>
        (
            IOobject::groupName(TName_, phaseName_)
        );

    tmp<volScalarField> mask =
        volScalarField::New("mask", mesh(), dimensionedScalar(dimless, 0));
    UIndirectList<scalar>(mask.ref().primitiveFieldRef(), zone_.zone()) = 1;
    const volScalarField htcAv
    (
        alpha*mask*heatTransferCoefficientModel_->htc()*heatTransferAv_->Av()
    );

    if (semiImplicit_)
    {
        if (he.dimensions() == dimEnergy/dimMass)
        {
            const basicThermo& thermo =
               mesh().lookupObject<basicThermo>
               (
                   IOobject::groupName(physicalProperties::typeName, phaseName_)
               );

            const volScalarField htcAvByCpv(htcAv/thermo.Cpv());

            eqn += htcAv*(Ta_ - T) + htcAvByCpv*he - fvm::Sp(htcAvByCpv, he);
        }
        else if (he.dimensions() == dimTemperature)
        {
            eqn += htcAv*Ta_ - fvm::Sp(htcAv, he);
        }
    }
    else
    {
        eqn += htcAv*(Ta_ - T);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransfer::heatTransfer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    zone_(mesh, dict),
    phaseName_(word::null),
    semiImplicit_(false),
    TName_(word::null),
    Ta_("Ta", dimTemperature, NaN),
    heatTransferAv_(nullptr),
    heatTransferCoefficientModel_(nullptr)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransfer::~heatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::heatTransfer::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    return wordList(1, thermo.he().name());
}


void Foam::fv::heatTransfer::addSup
(
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    add(geometricOneField(), eqn);
}


void Foam::fv::heatTransfer::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    add(geometricOneField(), eqn);
}


void Foam::fv::heatTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    add(alpha, eqn);
}


void Foam::fv::heatTransfer::correct()
{
    heatTransferCoefficientModel_->correct();
}


bool Foam::fv::heatTransfer::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::heatTransfer::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
}


void Foam::fv::heatTransfer::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::heatTransfer::distribute(const polyDistributionMap& map)
{
    zone_.distribute(map);
}


bool Foam::fv::heatTransfer::read(const dictionary& dict)
{
    if (fvModel::read(dict))
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

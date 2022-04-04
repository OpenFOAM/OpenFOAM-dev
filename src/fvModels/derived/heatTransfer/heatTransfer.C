/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

void Foam::fv::heatTransfer::readCoeffs()
{
    semiImplicit_ = coeffs().lookup<bool>("semiImplicit");

    TName_ = coeffs().lookupOrDefault<word>("T", "T");

    Ta_ = dimensionedScalar("Ta", dimTemperature, coeffs());

    heatTransferModel_ = heatTransferModel::New(coeffs(), mesh());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransfer::heatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    semiImplicit_(false),
    TName_(word::null),
    Ta_("Ta", dimTemperature, NaN),
    heatTransferModel_(nullptr)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransfer::~heatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::heatTransfer::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::heatTransfer::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const volScalarField& he = eqn.psi();

    const volScalarField& T =
        mesh().lookupObject<volScalarField>(TName_);

    tmp<volScalarField> mask =
        volScalarField::New("mask", mesh(), dimensionedScalar(dimless, 0));
    UIndirectList<scalar>(mask.ref().primitiveFieldRef(), set_.cells()) = 1;
    const volScalarField htcAoV
    (
        mask*heatTransferModel_->htc()*heatTransferModel_->AoV()
    );

    if (semiImplicit_)
    {
        if (he.dimensions() == dimEnergy/dimMass)
        {
            const basicThermo& thermo =
               mesh().lookupObject<basicThermo>(physicalProperties::typeName);

            const volScalarField htcAoVByCpv(htcAoV/thermo.Cpv());

            eqn += htcAoV*(Ta_ - T) + htcAoVByCpv*he - fvm::Sp(htcAoVByCpv, he);
        }
        else if (he.dimensions() == dimTemperature)
        {
            eqn += htcAoV*Ta_ - fvm::Sp(htcAoV, he);
        }
    }
    else
    {
        eqn += htcAoV*(Ta_ - T);
    }
}


void Foam::fv::heatTransfer::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    addSup(eqn, fieldName);
}


void Foam::fv::heatTransfer::correct()
{
    heatTransferModel_->correct();
}


bool Foam::fv::heatTransfer::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::heatTransfer::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::heatTransfer::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::heatTransfer::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::heatTransfer::read(const dictionary& dict)
{
    if (fvModel::read(dict))
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

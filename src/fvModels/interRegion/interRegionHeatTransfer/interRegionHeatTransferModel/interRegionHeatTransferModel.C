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

#include "interRegionHeatTransferModel.H"
#include "basicThermo.H"
#include "fvmSup.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcVolumeIntegrate.H"
#include "fvModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionHeatTransferModel, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::interRegionHeatTransferModel::readCoeffs()
{
    nbrModelName_ = coeffs().lookup<word>("nbrModel");

    semiImplicit_ = coeffs().lookup<bool>("semiImplicit");

    TName_ = coeffs().lookupOrDefault<word>("T", "T");
    TNbrName_ = coeffs().lookupOrDefault<word>("TNbr", "T");
}


Foam::fv::interRegionHeatTransferModel&
Foam::fv::interRegionHeatTransferModel::nbrModel() const
{
    const fvMesh& nbrMesh = mesh().time().lookupObject<fvMesh>(nbrRegionName());

    const PtrListDictionary<fvModel>& fvModels =
        nbrMesh.lookupObject<Foam::fvModels>("fvModels");

    if (fvModels.found(nbrModelName_))
    {
        return const_cast<interRegionHeatTransferModel&>
        (
            refCast<const interRegionHeatTransferModel>
            (
                fvModels[nbrModelName_]
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Neighbour model not found" << nbrModelName_
            << " in region " << nbrMesh.name() << nl
            << exit(FatalError);

        return const_cast<interRegionHeatTransferModel&>
        (
            NullObjectRef<interRegionHeatTransferModel>()
        );
    }
}


void Foam::fv::interRegionHeatTransferModel::correct() const
{
    if (master())
    {
        if (mesh().time().timeIndex() != timeIndex_)
        {
            correctHtc();
            timeIndex_ = mesh().time().timeIndex();
        }
    }
    else
    {
        nbrModel().correctHtc();
        interpolate(nbrModel().htc(), htc_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransferModel::interRegionHeatTransferModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionModel
    (
        name,
        modelType,
        dict,
        mesh
    ),
    nbrModelName_(word::null),
    timeIndex_(-1),
    semiImplicit_(false),
    TName_(word::null),
    TNbrName_(word::null),
    htc_
    (
        IOobject
        (
            type() + ":htc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimEnergy/dimTime/dimTemperature/dimVolume,
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    )
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransferModel::~interRegionHeatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::interRegionHeatTransferModel::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(basicThermo::dictName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::interRegionHeatTransferModel::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    correct();

    const volScalarField& he = eqn.psi();

    const volScalarField& T = mesh().lookupObject<volScalarField>(TName_);

    tmp<volScalarField> tTmapped
    (
        volScalarField::New(type() + ":Tmapped", T)
    );

    volScalarField& Tmapped = tTmapped.ref();

    const fvMesh& nbrMesh = mesh().time().lookupObject<fvMesh>(nbrRegionName());

    const volScalarField& Tnbr =
        nbrMesh.lookupObject<volScalarField>(TNbrName_);

    interpolate(Tnbr, Tmapped.primitiveFieldRef());

    if (debug)
    {
        Info<< "Volumetric integral of htc: "
            << fvc::domainIntegrate(htc_).value()
            << endl;

        if (mesh().time().writeTime())
        {
            Tmapped.write();
            htc_.write();
        }
    }

    if (semiImplicit_)
    {
        if (he.dimensions() == dimEnergy/dimMass)
        {
            const basicThermo& thermo =
               mesh().lookupObject<basicThermo>(basicThermo::dictName);

            const volScalarField htcByCpv(htc_/thermo.Cpv());

            eqn += htc_*(Tmapped - T) + htcByCpv*he - fvm::Sp(htcByCpv, he);

            if (debug)
            {
                const dimensionedScalar energy =
                    fvc::domainIntegrate(htc_*(Tmapped - T));

                Info<< "Energy exchange from region " << nbrMesh.name()
                    << " To " << mesh().name() << " : " <<  energy.value()
                    << endl;
            }
        }
        else if (he.dimensions() == dimTemperature)
        {
            eqn += htc_*Tmapped - fvm::Sp(htc_, he);
        }
    }
    else
    {
        eqn += htc_*(Tmapped - T);
    }
}


void Foam::fv::interRegionHeatTransferModel::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    addSup(eqn, fieldName);
}


bool Foam::fv::interRegionHeatTransferModel::read(const dictionary& dict)
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

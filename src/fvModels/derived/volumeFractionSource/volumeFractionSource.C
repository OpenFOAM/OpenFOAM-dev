/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "volumeFractionSource.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fluidThermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(volumeFractionSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        volumeFractionSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::volumeFractionSource::readCoeffs()
{
    phiName_ = coeffs().lookupOrDefault<word>("phi", "phi");
    rhoName_ = coeffs().lookupOrDefault<word>("rho", "rho");
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    volumePhaseName_ = coeffs().lookup<word>("volumePhase");
}


const Foam::volScalarField& Foam::fv::volumeFractionSource::volumeAlpha() const
{
    const word alphaName = IOobject::groupName("alpha", volumePhaseName_);

    if (!mesh().foundObject<volScalarField>(alphaName))
    {
        volScalarField* alphaPtr =
            new volScalarField
            (
                IOobject
                (
                    alphaName,
                    mesh().time().constant(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh()
            );

        alphaPtr->store();
    }

    return mesh().lookupObject<volScalarField>(alphaName);
}


Foam::tmp<Foam::volScalarField> Foam::fv::volumeFractionSource::D
(
    const word& fieldName
) const
{
    const word phiName =
        IOobject::groupName(phiName_, IOobject::group(fieldName));
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName);

    if (phi.dimensions() == dimVolume/dimTime)
    {
        const momentumTransportModel& turbulence =
            mesh().lookupType<momentumTransportModel>();

        return turbulence.nuEff();
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const fluidThermophysicalTransportModel& ttm =
            mesh().lookupType<fluidThermophysicalTransportModel>();

        return
            fieldName == ttm.thermo().T().name()
          ? ttm.kappaEff()
          : fieldName == ttm.thermo().he().name()
          ? ttm.kappaEff()/ttm.thermo().Cpv()
          : ttm.momentumTransport().rho()*ttm.momentumTransport().nuEff();
    }
    else
    {
        FatalErrorInFunction
            << "Dimensions of " << phi.name() << " not recognised"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


template <class Type, class AlphaFieldType>
void Foam::fv::volumeFractionSource::addGeneralSup
(
    const AlphaFieldType& alpha,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const word phiName =
        IOobject::groupName(phiName_, IOobject::group(fieldName));
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName);

    const volScalarField B(1 - volumeAlpha());
    const volScalarField AByB(volumeAlpha()/B);
    const volScalarField D(this->D(fieldName));

    // Divergence term
    const word divScheme = "div(" + phiName + "," + eqn.psi().name() + ")";
    eqn -= AByB*fvm::div(phi, eqn.psi(), divScheme);

    // Laplacian term
    const word laplacianScheme =
        "laplacian(" + D.name() + "," + eqn.psi().name() + ")";
    eqn +=
        fvm::laplacian(D, eqn.psi())
      - 1/B*fvm::laplacian(B*D, eqn.psi(), laplacianScheme);
}


template<class Type, class AlphaFieldType>
void Foam::fv::volumeFractionSource::addAlphaSupType
(
    const AlphaFieldType& alpha,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addGeneralSup(alpha, eqn, fieldName);
}


template<class AlphaFieldType>
void Foam::fv::volumeFractionSource::addAlphaSupType
(
    const AlphaFieldType& alpha,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == rhoName_)
    {
        const word phiName =
            IOobject::groupName(phiName_, IOobject::group(fieldName));
        const surfaceScalarField& phi =
            mesh().lookupObject<surfaceScalarField>(phiName);

        const volScalarField AByB(volumeAlpha()/(1 - volumeAlpha()));

        eqn -= AByB*fvc::div(phi);
    }
    else
    {
        addGeneralSup(alpha, eqn, fieldName);
    }
}


template<class AlphaFieldType>
void Foam::fv::volumeFractionSource::addAlphaSupType
(
    const AlphaFieldType& alpha,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == UName_)
    {
        const word phiName =
            IOobject::groupName(phiName_, IOobject::group(fieldName));
        const surfaceScalarField& phi =
            mesh().lookupObject<surfaceScalarField>(phiName);

        const volScalarField AByB(volumeAlpha()/(1 - volumeAlpha()));

        const word scheme("div(" + phiName + "," + eqn.psi().name() + ")");

        eqn -= fvm::div(fvc::interpolate(AByB)*phi, eqn.psi(), scheme);
    }
    else
    {
        addGeneralSup(alpha, eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::volumeFractionSource::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addAlphaSupType(geometricOneField(), eqn, fieldName);
}


template<class Type>
void Foam::fv::volumeFractionSource::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addAlphaSupType(geometricOneField(), eqn, fieldName);
}


template<class Type>
void Foam::fv::volumeFractionSource::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addAlphaSupType(alpha, eqn, fieldName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumeFractionSource::volumeFractionSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    phiName_(word::null),
    rhoName_(word::null),
    UName_(word::null),
    volumePhaseName_(word::null)
{
    readCoeffs();
    volumeAlpha();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::volumeFractionSource::~volumeFractionSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::volumeFractionSource::addsSupToField(const word& fieldName) const
{
    return true;
}


Foam::wordList Foam::fv::volumeFractionSource::addSupFields() const
{
    return wordList();
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_SUP, fv::volumeFractionSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_SUP, fv::volumeFractionSource);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_SUP,
    fv::volumeFractionSource
);


bool Foam::fv::volumeFractionSource::movePoints()
{
    return true;
}


void Foam::fv::volumeFractionSource::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::volumeFractionSource::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::volumeFractionSource::distribute(const polyDistributionMap&)
{}


bool Foam::fv::volumeFractionSource::read(const dictionary& dict)
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

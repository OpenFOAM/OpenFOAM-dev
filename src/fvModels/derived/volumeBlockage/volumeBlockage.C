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

#include "volumeBlockage.H"
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
    defineTypeNameAndDebug(volumeBlockage, 0);
    addToRunTimeSelectionTable(fvModel, volumeBlockage, dictionary);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvModel,
        volumeBlockage,
        dictionary,
        volumeFractionSource,
        "volumeFractionSource"
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::volumeBlockage::readCoeffs()
{
    phiName_ = coeffs().lookupOrDefault<word>("phi", "phi");
    rhoName_ = coeffs().lookupOrDefault<word>("rho", "rho");
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    volumePhaseName_ = coeffs().lookup<word>("volumePhase");
}


const Foam::volScalarField& Foam::fv::volumeBlockage::volumeAlpha() const
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


Foam::tmp<Foam::volScalarField> Foam::fv::volumeBlockage::D
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
void Foam::fv::volumeBlockage::addGeneralSupType
(
    const AlphaFieldType& alpha,
    fvMatrix<Type>& eqn
) const
{
    const word phiName =
        IOobject::groupName(phiName_, IOobject::group(eqn.psi().name()));
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName);

    const volScalarField B(1 - volumeAlpha());
    const volScalarField AByB(volumeAlpha()/B);
    const volScalarField D(this->D(eqn.psi().name()));

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
void Foam::fv::volumeBlockage::addAlphaSupType
(
    const AlphaFieldType& alpha,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    addGeneralSupType(alpha, eqn);
}


template<class AlphaFieldType>
void Foam::fv::volumeBlockage::addAlphaSupType
(
    const AlphaFieldType& alpha,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (IOobject::member(field.name()) == rhoName_)
    {
        const word phiName =
            IOobject::groupName(phiName_, IOobject::group(field.name()));
        const surfaceScalarField& phi =
            mesh().lookupObject<surfaceScalarField>(phiName);

        const volScalarField AByB(volumeAlpha()/(1 - volumeAlpha()));

        eqn -= AByB*fvc::div(phi);
    }
    else
    {
        addGeneralSupType(alpha, eqn);
    }
}


template<class AlphaFieldType>
void Foam::fv::volumeBlockage::addAlphaSupType
(
    const AlphaFieldType& alpha,
    const volVectorField& field,
    fvMatrix<vector>& eqn
) const
{
    if (IOobject::member(field.name()) == UName_)
    {
        const word phiName =
            IOobject::groupName(phiName_, IOobject::group(field.name()));
        const surfaceScalarField& phi =
            mesh().lookupObject<surfaceScalarField>(phiName);

        const volScalarField AByB(volumeAlpha()/(1 - volumeAlpha()));

        const word scheme("div(" + phiName + "," + eqn.psi().name() + ")");

        eqn -= fvm::div(fvc::interpolate(AByB)*phi, eqn.psi(), scheme);
    }
    else
    {
        addGeneralSupType(alpha, eqn);
    }
}


template<class Type>
void Foam::fv::volumeBlockage::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    addAlphaSupType(geometricOneField(), field, eqn);
}


template<class Type>
void Foam::fv::volumeBlockage::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    addAlphaSupType(geometricOneField(), field, eqn);
}


template<class Type>
void Foam::fv::volumeBlockage::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    addAlphaSupType(alpha, field, eqn);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumeBlockage::volumeBlockage
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

Foam::fv::volumeBlockage::~volumeBlockage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::volumeBlockage::addsSupToField(const word& fieldName) const
{
    return true;
}


Foam::wordList Foam::fv::volumeBlockage::addSupFields() const
{
    return wordList();
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_FIELD_SUP, fv::volumeBlockage)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP, fv::volumeBlockage)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::volumeBlockage
)


bool Foam::fv::volumeBlockage::movePoints()
{
    return true;
}


void Foam::fv::volumeBlockage::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::volumeBlockage::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::volumeBlockage::distribute(const polyDistributionMap&)
{}


bool Foam::fv::volumeBlockage::read(const dictionary& dict)
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

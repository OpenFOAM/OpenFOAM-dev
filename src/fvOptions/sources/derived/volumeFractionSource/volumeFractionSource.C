/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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
#include "thermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(volumeFractionSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        volumeFractionSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::volumeFractionSource::readCoeffs()
{
    phiName_ = coeffs_.lookupOrDefault<word>("phi", "phi");
    rhoName_ = coeffs_.lookupOrDefault<word>("rho", "rho");
    UName_ = coeffs_.lookupOrDefault<word>("U", "U");

    volumePhaseName_ = coeffs_.lookup<word>("volumePhase");
}


const Foam::volScalarField& Foam::fv::volumeFractionSource::volumeAlpha() const
{
    const word alphaName = IOobject::groupName("alpha", volumePhaseName_);

    if (!mesh_.foundObject<volScalarField>(alphaName))
    {
        volScalarField* alphaPtr =
            new volScalarField
            (
                IOobject
                (
                    alphaName,
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

        alphaPtr->store();
    }

    return mesh_.lookupObject<volScalarField>(alphaName);
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
            mesh().lookupObject<momentumTransportModel>
            (
                momentumTransportModel::typeName
            );

        return turbulence.nuEff();
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const thermophysicalTransportModel& ttm =
            mesh().lookupObject<thermophysicalTransportModel>
            (
                thermophysicalTransportModel::typeName
            );

        return
            fieldName == ttm.thermo().T().name()
          ? ttm.kappaEff()
          : fieldName == ttm.thermo().he().name()
          ? ttm.alphaEff()
          : ttm.momentumTransport().muEff();
    }
    else
    {
        FatalErrorInFunction
            << "Dimensions of " << phi.name() << " not recognised"
            << exit(FatalError);
        return tmp<volScalarField>(nullptr);
    }
}


template <class Type>
void Foam::fv::volumeFractionSource::addDivSup
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const word phiName =
        IOobject::groupName(phiName_, IOobject::group(fieldName));
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName);

    const volScalarField AByB(volumeAlpha()/(1 - volumeAlpha()));

    const word scheme("div(" + phiName + "," + eqn.psi().name() + ")");

    eqn -= AByB*fvm::div(phi, eqn.psi(), scheme);
}


void Foam::fv::volumeFractionSource::addUDivSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const word phiName =
        IOobject::groupName(phiName_, IOobject::group(fieldName));
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName);

    const volScalarField AByB(volumeAlpha()/(1 - volumeAlpha()));

    const word scheme("div(" + phiName + "," + eqn.psi().name() + ")");

    eqn -= fvm::div(fvc::interpolate(AByB)*phi, eqn.psi(), scheme);
}


void Foam::fv::volumeFractionSource::addRhoDivSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const word phiName =
        IOobject::groupName(phiName_, IOobject::group(fieldName));
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName);

    const volScalarField AByB(volumeAlpha()/(1 - volumeAlpha()));

    eqn -= AByB*fvc::div(phi);
}


template <class Type, class AlphaFieldType>
void Foam::fv::volumeFractionSource::addLaplacianSup
(
    const AlphaFieldType& alpha,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const volScalarField B(1 - volumeAlpha());

    const volScalarField D(this->D(fieldName));

    const word scheme("laplacian(" + D.name() + "," + eqn.psi().name() + ")");

    eqn +=
        fvm::laplacian(D, eqn.psi())
      - 1/B*fvm::laplacian(B*D, eqn.psi(), scheme);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumeFractionSource::volumeFractionSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
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

bool Foam::fv::volumeFractionSource::addsToField(const word& fieldName) const
{
    return true;
}


Foam::wordList Foam::fv::volumeFractionSource::addedToFields() const
{
    return wordList();
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == rhoName_)
    {
        addRhoDivSup(eqn, fieldName);
    }
    else
    {
        addDivSup(eqn, fieldName);
        addLaplacianSup(geometricOneField(), eqn, fieldName);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == UName_)
    {
        addUDivSup(eqn, fieldName);
    }
    else
    {
        addDivSup(eqn, fieldName);
        addLaplacianSup(geometricOneField(), eqn, fieldName);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(geometricOneField(), eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<symmTensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(geometricOneField(), eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<tensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(geometricOneField(), eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == rhoName_)
    {
        addRhoDivSup(eqn, fieldName);
    }
    else
    {
        addDivSup(eqn, fieldName);
        addLaplacianSup(geometricOneField(), eqn, fieldName);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == UName_)
    {
        addUDivSup(eqn, fieldName);
    }
    else
    {
        addDivSup(eqn, fieldName);
        addLaplacianSup(geometricOneField(), eqn, fieldName);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(geometricOneField(), eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(geometricOneField(), eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(geometricOneField(), eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == rhoName_)
    {
        addRhoDivSup(eqn, fieldName);
    }
    else
    {
        addDivSup(eqn, fieldName);
        addLaplacianSup(alpha, eqn, fieldName);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == UName_)
    {
        addUDivSup(eqn, fieldName);
    }
    else
    {
        addDivSup(eqn, fieldName);
        addLaplacianSup(alpha, eqn, fieldName);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(alpha, eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(alpha, eqn, fieldName);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const word& fieldName
) const
{
    addDivSup(eqn, fieldName);
    addLaplacianSup(alpha, eqn, fieldName);
}


bool Foam::fv::volumeFractionSource::read(const dictionary& dict)
{
    if (option::read(dict))
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

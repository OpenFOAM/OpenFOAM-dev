/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

const Foam::volScalarField& Foam::fv::volumeFractionSource::alpha() const
{
    const word alphaName = IOobject::groupName("alpha", phaseName_);

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
    const label fieldi
) const
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

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
            fieldNames_[fieldi] == ttm.thermo().T().name()
          ? ttm.kappaEff()
          : fieldNames_[fieldi] == ttm.thermo().he().name()
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
    const label fieldi
) const
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    const volScalarField AByB(this->alpha()/(1 - this->alpha()));

    eqn -= AByB*fvm::div(phi, eqn.psi());
}


void Foam::fv::volumeFractionSource::addUDivSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    const volScalarField AByB(this->alpha()/(1 - this->alpha()));

    const word scheme("div(" + phiName_ + "," + eqn.psi().name() + ")");

    eqn -= fvm::div(fvc::interpolate(AByB)*phi, eqn.psi(), scheme);
}


void Foam::fv::volumeFractionSource::addRhoDivSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    const volScalarField AByB(this->alpha()/(1 - this->alpha()));

    eqn -= AByB*fvc::div(phi);
}


template <class Type, class AlphaFieldType>
void Foam::fv::volumeFractionSource::addLaplacianSup
(
    const AlphaFieldType& alpha,
    fvMatrix<Type>& eqn,
    const label fieldi
) const
{
    const volScalarField B(1 - this->alpha());

    const volScalarField D(this->D(fieldi));

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
    phaseName_(dict.lookup<word>("phase")),
    phiName_("phi"),
    rhoName_("rho"),
    UName_("U")
{
    read(dict);
    alpha();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::volumeFractionSource::~volumeFractionSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    if (fieldNames_[fieldi] == rhoName_)
    {
        addRhoDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    if (fieldNames_[fieldi] == UName_)
    {
        addUDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    if (fieldNames_[fieldi] == rhoName_)
    {
        addRhoDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    if (fieldNames_[fieldi] == UName_)
    {
        addUDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    if (fieldNames_[fieldi] == rhoName_)
    {
        addRhoDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(alpha, eqn, fieldi);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    if (fieldNames_[fieldi] == UName_)
    {
        addUDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(alpha, eqn, fieldi);
    }
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(alpha, eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(alpha, eqn, fieldi);
}


void Foam::fv::volumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
) const
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(alpha, eqn, fieldi);
}


bool Foam::fv::volumeFractionSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        if (coeffs_.found("fields"))
        {
            coeffs_.lookup("fields") >> fieldNames_;
        }

        applied_.setSize(fieldNames_.size(), false);

        dict.readIfPresent("phi", phiName_);

        dict.readIfPresent("rho", rhoName_);

        dict.readIfPresent("U", UName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

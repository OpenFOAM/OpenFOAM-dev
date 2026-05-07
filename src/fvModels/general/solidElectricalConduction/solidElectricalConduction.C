/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "solidElectricalConduction.H"
#include "FunctionalGeometricField.H"
#include "basicThermo.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(solidElectricalConduction, 0);
    addToRunTimeSelectionTable(fvModel, solidElectricalConduction, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::solidElectricalConduction::readCoeffs(const dictionary& dict)
{
    const word sigmaKey = "sigma";
    const word sigmaScalarKey =
        sigmaKey + '<' + pTraits<scalar>::typeName + '>';
    const word sigmaTensorKey =
        sigmaKey + '<' + pTraits<tensor>::typeName + '>';

    const bool haveSigma = dict.found(sigmaKey);
    const bool haveScalarSigma = dict.found(sigmaScalarKey);
    const bool haveTensorSigma = dict.found(sigmaTensorKey);

    const label nHaveSigmas =
        label(haveSigma) + label(haveScalarSigma) + label(haveTensorSigma);

    if (nHaveSigmas != 1)
    {
        FatalIOErrorInFunction(dict)
            << (nHaveSigmas ? "multiple" : "none")
            << " of keywords " << sigmaKey << ", " << sigmaScalarKey << ' '
            << (nHaveSigmas ? "and" : "or") << ' ' << sigmaTensorKey
            << " defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }

    if (haveSigma || haveScalarSigma)
    {
        sigmaScalarPtr_.reset
        (
            new FunctionalGeometricField<scalar, fvMesh>
            (
                sigmaKey,
                haveSigma ? sigmaKey : sigmaScalarKey,
                mesh(),
                sqr(dimCurrent)/dimLength/dimPower,
                dict
            )
        );
    }
    else
    {
        sigmaTensorPtr_.reset
        (
            new FunctionalGeometricField<tensor, fvMesh>
            (
                sigmaKey,
                sigmaTensorKey,
                mesh(),
                sqr(dimCurrent)/dimLength/dimPower,
                dict
            )
        );
    }

    writeSigma_ = dict.lookupOrDefault<bool>("writeSigma", false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidElectricalConduction::solidElectricalConduction
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    phi_
    (
        IOobject
        (
            "phi",
            mesh().time().name(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    I_
    (
        IOobject
        (
            "I",
            mesh().time().name(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimCurrent, scalar(0))
    ),
    sigmaScalarPtr_(),
    sigmaTensorPtr_(),
    writeSigma_(false),
    meshChanged_(true)
{
    readCoeffs(dict);

    mesh().schemes().setFluxRequired(phi_.name());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::solidElectricalConduction::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::solidElectricalConduction::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    eqn += fvc::div(I_*fvc::interpolate(phi_));
}


void Foam::fv::solidElectricalConduction::correct()
{
    const bool sigmaChanged =
        sigmaScalarPtr_.valid()
      ? sigmaScalarPtr_->update()
      : sigmaTensorPtr_->update();

    if (sigmaChanged || meshChanged_)
    {
        const label nCorr =
            mesh()
           .solution()
           .solverDict(phi_.name())
           .lookupOrDefault<label>("nCorr", 0);

        for (label i = 0; i <= nCorr; ++ i)
        {
            fvMatrix<scalar> phiEqn
            (
                sigmaScalarPtr_.valid()
              ? fvm::laplacian(sigmaScalarPtr_(), phi_)
              : fvm::laplacian(sigmaTensorPtr_(), phi_)
            );

            phiEqn.solve();

            if (i == nCorr) I_ = phiEqn.flux();
        }
    }

    if (solutionControl::finalIteration(mesh()))
    {
        meshChanged_ = false;
    }
}


bool Foam::fv::solidElectricalConduction::movePoints()
{
    meshChanged_ = true;

    return true;
}


void Foam::fv::solidElectricalConduction::topoChange(const polyTopoChangeMap&)
{
    sigmaScalarPtr_.valid()
  ? sigmaScalarPtr_->reset()
  : sigmaTensorPtr_->reset();

    meshChanged_ = true;
}


void Foam::fv::solidElectricalConduction::mapMesh(const polyMeshMap&)
{
    sigmaScalarPtr_.valid()
  ? sigmaScalarPtr_->reset()
  : sigmaTensorPtr_->reset();

    meshChanged_ = true;
}


void Foam::fv::solidElectricalConduction::distribute(const polyDistributionMap&)
{
    sigmaScalarPtr_.valid()
  ? sigmaScalarPtr_->reset()
  : sigmaTensorPtr_->reset();

    meshChanged_ = true;
}


bool Foam::fv::solidElectricalConduction::read(const dictionary& dict)
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


bool Foam::fv::solidElectricalConduction::write(const bool write) const
{
    if (write && writeSigma_)
    {
        sigmaScalarPtr_.valid()
      ? sigmaScalarPtr_->write()
      : sigmaTensorPtr_->write();
    }

    return write;
}


// ************************************************************************* //

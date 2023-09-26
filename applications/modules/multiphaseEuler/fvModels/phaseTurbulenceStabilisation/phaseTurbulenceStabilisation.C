/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "phaseTurbulenceStabilisation.H"
#include "phaseSystem.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(phaseTurbulenceStabilisation, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            phaseTurbulenceStabilisation,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseTurbulenceStabilisation::addAlphaRhoSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn,
    tmp<volScalarField>
    (phaseCompressible::momentumTransportModel::*psi)() const
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const fvMesh& mesh = this->mesh();

    const phaseSystem::phaseModelPartialList& movingPhases =
        phase_.fluid().movingPhases();

    volScalarField::Internal transferRate
    (
        volScalarField::Internal::New
        (
            "transferRate",
            mesh,
            dimensionedScalar(dimless/dimTime, 0)
        )
    );
    volScalarField::Internal psiTransferRate
    (
        volScalarField::Internal::New
        (
            "psiTransferRate",
            mesh,
            dimensionedScalar((turbulence_.*psi)()().dimensions()/dimTime, 0)
        )
    );

    forAll(movingPhases, phasei)
    {
        if (movingPhases[phasei] != phase_)
        {
            const phaseCompressible::momentumTransportModel& turbulence =
                mesh.lookupType<phaseCompressible::momentumTransportModel>
                (
                    phaseName_
                );

            if (notNull(turbulence))
            {
                const volScalarField::Internal phaseTransferRate
                (
                    movingPhases[phasei]
                   *min
                    (
                        turbulence.epsilon()/turbulence.k(),
                        1.0/phase_.time().deltaT()
                    )
                );

                transferRate += phaseTransferRate;
                psiTransferRate += phaseTransferRate*(turbulence.*psi)()();
            }
        }
    }

    const volScalarField::Internal transferCoeff
    (
        max(alphaInversion_ - alpha(), scalar(0))*rho()
    );

    eqn += transferCoeff*psiTransferRate;
    eqn -= fvm::Sp(transferCoeff*transferRate, eqn.psi());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseTurbulenceStabilisation::phaseTurbulenceStabilisation
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    phaseName_(dict.lookup("phase")),
    alphaInversion_("alphaInversion", dimless, dict),
    phase_
    (
        mesh.lookupObject<phaseModel>(IOobject::groupName("alpha", phaseName_))
    ),
    turbulence_
    (
        mesh.lookupType<phaseCompressible::momentumTransportModel>
        (
            phaseName_
        )
    )
{
    const word kName(IOobject::groupName("k", phaseName_));
    const word epsilonName(IOobject::groupName("epsilon", phaseName_));
    const word omegaName(IOobject::groupName("omega", phaseName_));

    if (mesh.foundObject<volScalarField>(kName))
    {
        fieldNames_.append(kName);
    }

    if (mesh.foundObject<volScalarField>(epsilonName))
    {
        fieldNames_.append(epsilonName);
    }

    if (mesh.foundObject<volScalarField>(omegaName))
    {
        fieldNames_.append(omegaName);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::phaseTurbulenceStabilisation::addSupFields() const
{
    return fieldNames_;
}


void Foam::fv::phaseTurbulenceStabilisation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (field.name() == IOobject::groupName("k", phaseName_))
    {
        addAlphaRhoSup
        (
            alpha,
            rho,
            field,
            eqn,
            &phaseCompressible::momentumTransportModel::k
        );
    }
    else if (field.name() == IOobject::groupName("epsilon", phaseName_))
    {
        addAlphaRhoSup
        (
            alpha,
            rho,
            field,
            eqn,
            &phaseCompressible::momentumTransportModel::epsilon
        );
    }
    else if (field.name() == IOobject::groupName("omega", phaseName_))
    {
        addAlphaRhoSup
        (
            alpha,
            rho,
            field,
            eqn,
            &phaseCompressible::momentumTransportModel::omega
        );
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << field.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::phaseTurbulenceStabilisation::topoChange
(
    const polyTopoChangeMap&
)
{}


void Foam::fv::phaseTurbulenceStabilisation::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::phaseTurbulenceStabilisation::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fv::phaseTurbulenceStabilisation::movePoints()
{
    return true;
}


// ************************************************************************* //

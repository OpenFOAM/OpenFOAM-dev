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

#include "interfaceTurbulenceDamping.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(interfaceTurbulenceDamping, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            interfaceTurbulenceDamping,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::interfaceTurbulenceDamping::interfaceFraction
(
    const volScalarField& alpha
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<volScalarField::Internal> tA
    (
        volScalarField::Internal::New
        (
            "A",
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField::Internal& A = tA.ref();

    const surfaceVectorField& Sf = mesh.Sf();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const surfaceScalarField alphaf(fvc::interpolate(alpha));

    const volVectorField gradAlpha(fvc::grad(alpha));
    const volVectorField::Internal n
    (
        gradAlpha()/(mag(gradAlpha()) + phase_.fluid().deltaN())
    );

    const scalarField& ialpha = alpha;
    const scalarField& ialphaf = alphaf;
    scalarField sumnSf(mesh.nCells(), 0);

    forAll(own, facei)
    {
        {
            const scalar nSf(mag(n[own[facei]] & Sf[facei]));
            A[own[facei]] += nSf*(ialphaf[facei] - ialpha[own[facei]]);
            sumnSf[own[facei]] += nSf;
        }
        {
            const scalar nSf(mag(n[nei[facei]] & Sf[facei]));
            A[nei[facei]] += nSf*(ialphaf[facei] - ialpha[nei[facei]]);
            sumnSf[nei[facei]] += nSf;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& own = mesh.boundary()[patchi].faceCells();
        const fvsPatchScalarField& palphaf = alphaf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            const scalar nSf(mag(n[own[facei]] & Sf[facei]));
            A[own[facei]] += nSf*(palphaf[facei] - ialpha[own[facei]]);
            sumnSf[own[facei]] += nSf;
        }
    }

    scalarField& a = A.field();
    forAll(a, i)
    {
        if (sumnSf[i] > small)
        {
            a[i] = 2*mag(a[i])/sumnSf[i];
        }
        else
        {
            a[i] = 0;
        }
    }

    return tA;
}


template<class RhoType>
void Foam::fv::interfaceTurbulenceDamping::addRhoSup
(
    const RhoType& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const phaseSystem::phaseModelPartialList& movingPhases =
        phase_.fluid().movingPhases();

    volScalarField::Internal aSqrnu
    (
        movingPhases[0]*sqr(movingPhases[0].thermo().nu()()())
    );

    for (label phasei=1; phasei<movingPhases.size(); phasei++)
    {
        aSqrnu +=
            movingPhases[phasei]*sqr(movingPhases[phasei].thermo().nu()()());
    }

    if (fieldName == "epsilon")
    {
        eqn += rho*interfaceFraction(phase_)*C2_*aSqrnu*turbulence_.k()()
            /pow4(delta_);
    }
    else if (fieldName == "omega")
    {
        eqn += rho*interfaceFraction(phase_)*beta_*aSqrnu
            /(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interfaceTurbulenceDamping::interfaceTurbulenceDamping
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    phaseName_(dict.lookup("phase")),
    delta_("delta", dimLength, dict),
    phase_
    (
        mesh.lookupObject<phaseModel>(IOobject::groupName("alpha", phaseName_))
    ),
    turbulence_
    (
        mesh.lookupType<phaseCompressible::momentumTransportModel>(phaseName_)
    ),
    C2_("C2", dimless, 0),
    betaStar_("betaStar", dimless, 0),
    beta_("beta", dimless, 0)
{
    const word epsilonName(IOobject::groupName("epsilon", phaseName_));
    const word omegaName(IOobject::groupName("omega", phaseName_));

    if (mesh.foundObject<volScalarField>(epsilonName))
    {
        fieldName_ = epsilonName;
        C2_.read(turbulence_.coeffDict());
    }
    else if (mesh.foundObject<volScalarField>(omegaName))
    {
        fieldName_ = omegaName;
        betaStar_.read(turbulence_.coeffDict());

        // Read beta for k-omega models or beta1 for k-omega SST
        if (turbulence_.coeffDict().found("beta"))
        {
            beta_.read(turbulence_.coeffDict());
        }
        else
        {
            beta_ =
                dimensionedScalar("beta1", dimless, turbulence_.coeffDict());
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Cannot find either " << epsilonName << " or " << omegaName
            << " field for fvModel " << typeName << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::interfaceTurbulenceDamping::addSupFields() const
{
    return wordList(1, fieldName_);
}


void Foam::fv::interfaceTurbulenceDamping::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    addRhoSup(one(), eqn, fieldName);
}


void Foam::fv::interfaceTurbulenceDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    addRhoSup(rho(), eqn, fieldName);
}


void Foam::fv::interfaceTurbulenceDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    volScalarField::Internal aSqrnu
    (
        alpha*sqr(phase_.thermo().nu()()())
    );

    if (fieldName == IOobject::groupName("epsilon", phaseName_))
    {
        eqn += rho()*interfaceFraction(alpha)
            *C2_*aSqrnu*turbulence_.k()()/pow4(delta_);
    }
    else if (fieldName == IOobject::groupName("omega", phaseName_))
    {
        eqn += rho()*interfaceFraction(alpha)
            *beta_*aSqrnu/(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::interfaceTurbulenceDamping::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::interfaceTurbulenceDamping::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::interfaceTurbulenceDamping::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fv::interfaceTurbulenceDamping::movePoints()
{
    return true;
}


// ************************************************************************* //

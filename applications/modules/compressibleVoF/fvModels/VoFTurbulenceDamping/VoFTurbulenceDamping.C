/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "VoFTurbulenceDamping.H"
#include "compressibleTwoPhaseVoFMixture.H"
#include "compressibleMomentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        namespace compressible
        {
            defineTypeNameAndDebug(VoFTurbulenceDamping, 0);

            addToRunTimeSelectionTable
            (
                fvModel,
                VoFTurbulenceDamping,
                dictionary
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::compressible::VoFTurbulenceDamping::VoFTurbulenceDamping
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    phaseName_(dict.lookupOrDefault("phase", word::null)),
    delta_("delta", dimLength, dict),
    mixture_
    (
        mesh.lookupObject<compressibleTwoPhaseVoFMixture>
        (
            "phaseProperties"
        )
    ),
    momentumTransport_
    (
        mesh.lookupType<compressibleMomentumTransportModel>(phaseName_)
    ),
    C2_("C2", dimless, dict, 1.92),
    betaStar_("betaStar", dimless, dict, 0.09),
    beta_("beta", dimless, dict, 0.072)
{
    const word epsilonName(IOobject::groupName("epsilon", phaseName_));
    const word omegaName(IOobject::groupName("omega", phaseName_));

    if (mesh.foundObject<volScalarField>(epsilonName))
    {
        fieldName_ = epsilonName;
    }
    else if (mesh.foundObject<volScalarField>(omegaName))
    {
        fieldName_ = omegaName;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Cannot find either " << epsilonName << " or " << omegaName
            << " field for fvModel " << typeName << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::fv::compressible::VoFTurbulenceDamping::addSupFields() const
{
    return wordList(1, fieldName_);
}


void Foam::fv::compressible::VoFTurbulenceDamping::addSup
(
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField::Internal aRhoSqrnu
    (
        mixture_.alpha1()()*mixture_.rho1()()*sqr(mixture_.thermo1().nu()()())
      + mixture_.alpha2()()*mixture_.rho2()()*sqr(mixture_.thermo2().nu()()())
    );

    if (field.name() == "epsilon")
    {
        eqn += mixture_.interfaceFraction()
           *C2_*aRhoSqrnu*momentumTransport_.k()()/pow4(delta_);
    }
    else if (field.name() == "omega")
    {
        eqn += mixture_.interfaceFraction()
           *beta_*aRhoSqrnu/(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << field.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::compressible::VoFTurbulenceDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    tmp<volScalarField::Internal> taRhoSqrnu;

    if (mixture_.alpha1().name() == alpha.name())
    {
        taRhoSqrnu = mixture_.alpha1()()*mixture_.rho1()()
           *sqr(mixture_.thermo1().nu()()());
    }
    else if (mixture_.alpha2().name() == alpha.name())
    {
        taRhoSqrnu = mixture_.alpha2()()*mixture_.rho2()()
           *sqr(mixture_.thermo2().nu()()());
    }
    else
    {
        FatalErrorInFunction
            << "Unknown phase-fraction " << alpha.name()
            << exit(FatalError);
    }

    if (field.name() == IOobject::groupName("epsilon", phaseName_))
    {
        eqn += mixture_.interfaceFraction()
           *C2_*taRhoSqrnu*momentumTransport_.k()()/pow4(delta_);
    }
    else if (field.name() == IOobject::groupName("omega", phaseName_))
    {
        eqn += mixture_.interfaceFraction()
           *beta_*taRhoSqrnu/(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << field.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::compressible::VoFTurbulenceDamping::topoChange
(
    const polyTopoChangeMap&
)
{}


void Foam::fv::compressible::VoFTurbulenceDamping::mapMesh
(
    const polyMeshMap& map
)
{}


void Foam::fv::compressible::VoFTurbulenceDamping::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fv::compressible::VoFTurbulenceDamping::movePoints()
{
    return true;
}


// ************************************************************************* //

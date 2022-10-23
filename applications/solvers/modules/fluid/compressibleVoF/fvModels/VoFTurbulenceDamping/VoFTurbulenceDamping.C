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

#include "VoFTurbulenceDamping.H"
#include "compressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"
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
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(sourceName, modelType, dict, mesh),
    phaseName_(dict.lookupOrDefault("phase", word::null)),
    delta_("delta", dimLength, dict),
    mixture_
    (
        mesh.lookupObject<compressibleTwoPhaseMixture>
        (
            "phaseProperties"
        )
    ),
    interface_(mixture_),
    turbulence_
    (
        mesh.lookupType<compressibleMomentumTransportModel>(phaseName_)
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

Foam::wordList
Foam::fv::compressible::VoFTurbulenceDamping::addSupFields() const
{
    return wordList(1, fieldName_);
}


void Foam::fv::compressible::VoFTurbulenceDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
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

    if (fieldName == "epsilon")
    {
        eqn += interface_.fraction()*C2_*aRhoSqrnu*turbulence_.k()()
           /pow4(delta_);
    }
    else if (fieldName == "omega")
    {
        eqn += interface_.fraction()*beta_*aRhoSqrnu
           /(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::compressible::VoFTurbulenceDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField&,
    fvMatrix<scalar>& eqn,
    const word& fieldName
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

    if (fieldName == IOobject::groupName("epsilon", phaseName_))
    {
        eqn += interface_.fraction()*C2_*taRhoSqrnu*turbulence_.k()()
           /pow4(delta_);
    }
    else if (fieldName == IOobject::groupName("omega", phaseName_))
    {
        eqn += interface_.fraction()*beta_*taRhoSqrnu
           /(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
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

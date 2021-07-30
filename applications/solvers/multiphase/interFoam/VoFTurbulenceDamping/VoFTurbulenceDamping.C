/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
#include "incompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"
#include "incompressibleMomentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFTurbulenceDamping::VoFTurbulenceDamping
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
        mesh.lookupObject<incompressibleTwoPhaseMixture>
        (
            "phaseProperties"
        )
    ),
    interface_(refCast<const interfaceProperties>(mixture_)),
    turbulence_
    (
        mesh.lookupObject<incompressibleMomentumTransportModel>
        (IOobject::groupName(momentumTransportModel::typeName, phaseName_))
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

Foam::wordList Foam::fv::VoFTurbulenceDamping::addSupFields() const
{
    return wordList(1, fieldName_);
}


void Foam::fv::VoFTurbulenceDamping::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField::Internal aSqrnu
    (
        mixture_.alpha1()()*sqr(mixture_.nuModel1().nu()()())
      + mixture_.alpha2()()*sqr(mixture_.nuModel2().nu()()())
    );

    if (fieldName == "epsilon")
    {
        eqn += interface_.fraction()*C2_*aSqrnu*turbulence_.k()()/pow4(delta_);
    }
    else if (fieldName == "omega")
    {
        eqn += interface_.fraction()*beta_*aSqrnu/(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFTurbulenceDamping::addSup
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

    tmp<volScalarField::Internal> taSqrnu;

    if (mixture_.alpha1().name() == alpha.name())
    {
        taSqrnu = mixture_.alpha1()()*sqr(mixture_.nuModel1().nu()()());
    }
    else if (mixture_.alpha2().name() == alpha.name())
    {
        taSqrnu = mixture_.alpha2()()*sqr(mixture_.nuModel2().nu()()());
    }
    else
    {
        FatalErrorInFunction
            << "Unknown phase-fraction " << alpha.name()
            << exit(FatalError);
    }

    if (fieldName == IOobject::groupName("epsilon", phaseName_))
    {
        eqn += interface_.fraction()*C2_*taSqrnu*turbulence_.k()()
           /pow4(delta_);
    }
    else if (fieldName == IOobject::groupName("omega", phaseName_))
    {
        eqn += interface_.fraction()*beta_*taSqrnu
           /(sqr(betaStar_)*pow4(delta_));
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


// ************************************************************************* //

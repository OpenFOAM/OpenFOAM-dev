/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "Gulder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(Gulder, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        Gulder,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::Gulder::Gulder
(
    const dictionary& dict,
    const dictionary& coeffDict,
    const psiuMulticomponentThermo& ct
)
:
    laminarFlameSpeed(dict, ct),

    W_(coeffDict.lookup<scalar>("W")),
    eta_(coeffDict.lookup<scalar>("eta")),
    xi_(coeffDict.lookup<scalar>("xi")),
    f_(coeffDict.lookup<scalar>("f")),
    alpha_(coeffDict.lookup<scalar>("alpha")),
    beta_(coeffDict.lookup<scalar>("beta"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::Gulder::~Gulder()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::laminarFlameSpeedModels::Gulder::SuRef
(
    const scalar Phi
) const
{
    if (Phi > small)
    {
        return W_*pow(Phi, eta_)*exp(-xi_*sqr(Phi - 1.075));
    }
    else
    {
        return 0;
    }
}


inline Foam::scalar Foam::laminarFlameSpeedModels::Gulder::Su0pTphi
(
    const scalar p,
    const scalar Tu,
    const scalar Phi,
    const scalar Yegr
) const
{
    static const scalar Tref = 300.0;
    static const scalar pRef = 1.013e5;

    return
        SuRef(Phi)*pow((Tu/Tref), alpha_)*pow((p/pRef), beta_)
       *max(1 - f_*Yegr, 0);
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::Gulder::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    const volScalarField& phi
) const
{
    tmp<volScalarField> tSu0
    (
        volScalarField::New
        (
            "Su0",
            p.mesh(),
            dimensionedScalar(dimVelocity, 0)
        )
    );

    volScalarField& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi[celli], 0);
    }

    volScalarField::Boundary& Su0Bf = Su0.boundaryFieldRef();

    forAll(Su0Bf, patchi)
    {
        forAll(Su0Bf[patchi], facei)
        {
            Su0Bf[patchi][facei] =
                Su0pTphi
                (
                    p.boundaryField()[patchi][facei],
                    Tu.boundaryField()[patchi][facei],
                    phi.boundaryField()[patchi][facei],
                    0
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::Gulder::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    const volScalarField& Phi,
    const volScalarField& egr
) const
{
    tmp<volScalarField> tSu0
    (
        volScalarField::New
        (
            "Su0",
            p.mesh(),
            dimensionedScalar(dimVelocity, 0)
        )
    );

    volScalarField& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], Phi[celli], egr[celli]);
    }

    volScalarField::Boundary& Su0Bf = Su0.boundaryFieldRef();

    forAll(Su0Bf, patchi)
    {
        forAll(Su0Bf[patchi], facei)
        {
            Su0Bf[patchi][facei] =
                Su0pTphi
                (
                    p.boundaryField()[patchi][facei],
                    Tu.boundaryField()[patchi][facei],
                    Phi.boundaryField()[patchi][facei],
                    egr.boundaryField()[patchi][facei]
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::Gulder::operator()() const
{
    if (psiuMulticomponentThermo_.containsSpecie("egr"))
    {
        return Su0pTphi
        (
            psiuMulticomponentThermo_.p(),
            psiuMulticomponentThermo_.Tu(),
            psiuMulticomponentThermo_.Phi(),
            psiuMulticomponentThermo_.Y("egr")
        );
    }
    else
    {
        return Su0pTphi
        (
            psiuMulticomponentThermo_.p(),
            psiuMulticomponentThermo_.Tu(),
            psiuMulticomponentThermo_.Phi()
        );
    }
}


// ************************************************************************* //

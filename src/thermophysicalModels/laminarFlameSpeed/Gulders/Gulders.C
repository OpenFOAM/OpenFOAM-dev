/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "Gulders.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(Gulders, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        Gulders,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::Gulders::Gulders
(
    const dictionary& dict,
    const psiuMulticomponentThermo& ct
)
:
    laminarFlameSpeed(dict, ct),

    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs").subDict(fuel_)),
    W_(coeffsDict_.lookup<scalar>("W")),
    eta_(coeffsDict_.lookup<scalar>("eta")),
    xi_(coeffsDict_.lookup<scalar>("xi")),
    f_(coeffsDict_.lookup<scalar>("f")),
    alpha_(coeffsDict_.lookup<scalar>("alpha")),
    beta_(coeffsDict_.lookup<scalar>("beta"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::Gulders::~Gulders()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::laminarFlameSpeedModels::Gulders::SuRef
(
    scalar phi
) const
{
    if (phi > small)
    {
        return W_*pow(phi, eta_)*exp(-xi_*sqr(phi - 1.075));
    }
    else
    {
        return 0.0;
    }
}


inline Foam::scalar Foam::laminarFlameSpeedModels::Gulders::Su0pTphi
(
    scalar p,
    scalar Tu,
    scalar phi,
    scalar Yres
) const
{
    static const scalar Tref = 300.0;
    static const scalar pRef = 1.013e5;

    return SuRef(phi)*pow((Tu/Tref), alpha_)*pow((p/pRef), beta_)*(1 - f_*Yres);
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::Gulders::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    scalar phi
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
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi, 0.0);
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
                    phi,
                    0.0
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::Gulders::Su0pTphi
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
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi[celli], 0.0);
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
                    0.0
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::Gulders::operator()() const
{
    if (psiuMulticomponentThermo_.composition().contains("ft"))
    {
        const volScalarField& ft =
            psiuMulticomponentThermo_.composition().Y("ft");

        return Su0pTphi
        (
            psiuMulticomponentThermo_.p(),
            psiuMulticomponentThermo_.Tu(),
            dimensionedScalar
            (
                "stoichiometricAirFuelMassRatio",
                dimless,
                psiuMulticomponentThermo_.properties()
            )*ft/max(1 - ft, small)
        );
    }
    else
    {
        return Su0pTphi
        (
            psiuMulticomponentThermo_.p(),
            psiuMulticomponentThermo_.Tu(),
            equivalenceRatio_
        );
    }
}


// ************************************************************************* //

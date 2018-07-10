/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "IFstream.H"
#include "SCOPELaminarFlameSpeed.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(SCOPE, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        SCOPE,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::SCOPE::polynomial::polynomial
(
    const dictionary& polyDict
)
:
    FixedList<scalar, 7>(polyDict.lookup("coefficients")),
    ll(readScalar(polyDict.lookup("lowerLimit"))),
    ul(readScalar(polyDict.lookup("upperLimit"))),
    llv(polyPhi(ll, *this)),
    ulv(polyPhi(ul, *this)),
    lu(0)
{}


Foam::laminarFlameSpeedModels::SCOPE::SCOPE
(
    const dictionary& dict,
    const psiuReactionThermo& ct
)
:
    laminarFlameSpeed(dict, ct),

    coeffsDict_
    (
        dictionary
        (
          IFstream
          (
              fileName
              (
                  dict.lookup("fuelFile")
              )
          )()
        ).optionalSubDict(typeName + "Coeffs")
    ),
    LFL_(readScalar(coeffsDict_.lookup("lowerFlamabilityLimit"))),
    UFL_(readScalar(coeffsDict_.lookup("upperFlamabilityLimit"))),
    SuPolyL_(coeffsDict_.subDict("lowerSuPolynomial")),
    SuPolyU_(coeffsDict_.subDict("upperSuPolynomial")),
    Texp_(readScalar(coeffsDict_.lookup("Texp"))),
    pexp_(readScalar(coeffsDict_.lookup("pexp"))),
    MaPolyL_(coeffsDict_.subDict("lowerMaPolynomial")),
    MaPolyU_(coeffsDict_.subDict("upperMaPolynomial"))
{
    SuPolyL_.ll = max(SuPolyL_.ll, LFL_) + small;
    SuPolyU_.ul = min(SuPolyU_.ul, UFL_) - small;

    SuPolyL_.lu = 0.5*(SuPolyL_.ul + SuPolyU_.ll);
    SuPolyU_.lu = SuPolyL_.lu - small;

    MaPolyL_.lu = 0.5*(MaPolyL_.ul + MaPolyU_.ll);
    MaPolyU_.lu = MaPolyL_.lu - small;

    if (debug)
    {
        Info<< "phi     Su  (T = Tref, p = pref)" << endl;
        label n = 200;
        for (int i=0; i<n; i++)
        {
            scalar phi = (2.0*i)/n;
            Info<< phi << token::TAB << SuRef(phi) << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::SCOPE::~SCOPE()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::polyPhi
(
    scalar phi,
    const polynomial& a
)
{
    scalar x = phi - 1.0;

    return
        a[0]
       *(
           scalar(1)
         + x*(a[1] + x*(a[2] + x*(a[3] + x*(a[4] + x*(a[5] + x*a[6])))))
        );
}


inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::SuRef
(
    scalar phi
) const
{
    if (phi < LFL_ || phi > UFL_)
    {
        // Return 0 beyond the flamibility limits
        return scalar(0);
    }
    else if (phi < SuPolyL_.ll)
    {
        // Use linear interpolation between the low end of the
        // lower polynomial and the lower flamibility limit
        return SuPolyL_.llv*(phi - LFL_)/(SuPolyL_.ll - LFL_);
    }
    else if (phi > SuPolyU_.ul)
    {
        // Use linear interpolation between the upper end of the
        // upper polynomial and the upper flamibility limit
        return SuPolyU_.ulv*(UFL_ - phi)/(UFL_ - SuPolyU_.ul);
    }
    else if (phi < SuPolyL_.lu)
    {
        // Evaluate the lower polynomial
        return polyPhi(phi, SuPolyL_);
    }
    else if (phi > SuPolyU_.lu)
    {
        // Evaluate the upper polynomial
        return polyPhi(phi, SuPolyU_);
    }
    else
    {
        FatalErrorInFunction
            << "phi = " << phi
            << " cannot be handled by SCOPE function with the "
               "given coefficients"
            << exit(FatalError);

        return scalar(0);
    }
}


inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::Ma
(
    scalar phi
) const
{
    if (phi < MaPolyL_.ll)
    {
        // Beyond the lower limit assume Ma is constant
        return MaPolyL_.llv;
    }
    else if (phi > MaPolyU_.ul)
    {
        // Beyond the upper limit assume Ma is constant
        return MaPolyU_.ulv;
    }
    else if (phi < SuPolyL_.lu)
    {
        // Evaluate the lower polynomial
        return polyPhi(phi, MaPolyL_);
    }
    else if (phi > SuPolyU_.lu)
    {
        // Evaluate the upper polynomial
        return polyPhi(phi, MaPolyU_);
    }
    else
    {
        FatalErrorInFunction
            << "phi = " << phi
            << " cannot be handled by SCOPE function with the "
               "given coefficients"
            << exit(FatalError);

        return scalar(0);
    }
}


inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::Su0pTphi
(
    scalar p,
    scalar Tu,
    scalar phi
) const
{
    static const scalar Tref = 300.0;
    static const scalar pRef = 1.013e5;

    return SuRef(phi)*pow((Tu/Tref), Texp_)*pow((p/pRef), pexp_);
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::SCOPE::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    scalar phi
) const
{
    tmp<volScalarField> tSu0
    (
        new volScalarField
        (
            IOobject
            (
                "Su0",
                p.time().timeName(),
                p.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("Su0", dimVelocity, 0.0)
        )
    );

    volScalarField& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi);
    }

    volScalarField::Boundary& Su0Bf = Su0.boundaryFieldRef();

    forAll(Su0Bf, patchi)
    {
        scalarField& Su0p = Su0Bf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tup = Tu.boundaryField()[patchi];

        forAll(Su0p, facei)
        {
            Su0p[facei] = Su0pTphi(pp[facei], Tup[facei], phi);
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::SCOPE::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    const volScalarField& phi
) const
{
    tmp<volScalarField> tSu0
    (
        new volScalarField
        (
            IOobject
            (
                "Su0",
                p.time().timeName(),
                p.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("Su0", dimVelocity, 0.0)
        )
    );

    volScalarField& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi[celli]);
    }

    volScalarField::Boundary& Su0Bf = Su0.boundaryFieldRef();

    forAll(Su0Bf, patchi)
    {
        scalarField& Su0p = Su0Bf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tup = Tu.boundaryField()[patchi];
        const scalarField& phip = phi.boundaryField()[patchi];

        forAll(Su0p, facei)
        {
            Su0p[facei] =
                Su0pTphi
                (
                    pp[facei],
                    Tup[facei],
                    phip[facei]
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::SCOPE::Ma
(
    const volScalarField& phi
) const
{
    tmp<volScalarField> tMa
    (
        new volScalarField
        (
            IOobject
            (
                "Ma",
                phi.time().timeName(),
                phi.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi.mesh(),
            dimensionedScalar("Ma", dimless, 0.0)
        )
    );

    volScalarField& ma = tMa.ref();

    forAll(ma, celli)
    {
        ma[celli] = Ma(phi[celli]);
    }

    volScalarField::Boundary& maBf = ma.boundaryFieldRef();

    forAll(maBf, patchi)
    {
        scalarField& map = maBf[patchi];
        const scalarField& phip = phi.boundaryField()[patchi];

        forAll(map, facei)
        {
            map[facei] = Ma(phip[facei]);
        }
    }

    return tMa;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::SCOPE::Ma() const
{
    if (psiuReactionThermo_.composition().contains("ft"))
    {
        const volScalarField& ft = psiuReactionThermo_.composition().Y("ft");

        return Ma
        (
            dimensionedScalar
            (
                psiuReactionThermo_.lookup("stoichiometricAirFuelMassRatio")
            )*ft/(scalar(1) - ft)
        );
    }
    else
    {
        const fvMesh& mesh = psiuReactionThermo_.p().mesh();

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Ma",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("Ma", dimless, Ma(equivalenceRatio_))
            )
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::SCOPE::operator()() const
{
    if (psiuReactionThermo_.composition().contains("ft"))
    {
        const volScalarField& ft = psiuReactionThermo_.composition().Y("ft");

        return Su0pTphi
        (
            psiuReactionThermo_.p(),
            psiuReactionThermo_.Tu(),
            dimensionedScalar
            (
                psiuReactionThermo_.lookup("stoichiometricAirFuelMassRatio")
            )*ft/(scalar(1) - ft)
        );
    }
    else
    {
        return Su0pTphi
        (
            psiuReactionThermo_.p(),
            psiuReactionThermo_.Tu(),
            equivalenceRatio_
        );
    }
}


// ************************************************************************* //

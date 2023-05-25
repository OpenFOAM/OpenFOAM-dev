/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "segregated.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(segregated, 0);
    addToRunTimeSelectionTable(dragModel, segregated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::segregated::segregated
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dragModel(dict, interface, registerObject),
    interface_(interface.modelCast<dragModel, segregatedPhaseInterface>()),
    m_("m", dimless, dict),
    n_("n", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::segregated::~segregated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::segregated::K() const
{
    const fvMesh& mesh(interface_.phase1().mesh());

    const volScalarField& alpha1(interface_.phase1());
    const volScalarField& alpha2(interface_.phase2());

    const volScalarField::Internal& rho1(interface_.phase1().rho());
    const volScalarField::Internal& rho2(interface_.phase2().rho());

    tmp<volScalarField> tnu1(interface_.phase1().thermo().nu());
    tmp<volScalarField> tnu2(interface_.phase2().thermo().nu());

    const volScalarField::Internal& nu1(tnu1());
    const volScalarField::Internal& nu2(tnu2());

    const volScalarField::Internal L(cbrt(mesh.V()));

    const dimensionedScalar residualAlpha
    (
        (
            interface_.phase1().residualAlpha()
          + interface_.phase2().residualAlpha()
        )/2
    );

    const volScalarField I1
    (
        alpha1/max(alpha1 + alpha2, residualAlpha)
    );
    const volScalarField I2
    (
        alpha2/max(alpha1 + alpha2, residualAlpha)
    );
    const volScalarField::Internal magGradI
    (
        max
        (
            (
                rho2*mag(fvc::grad(I1)()())
              + rho1*mag(fvc::grad(I2)()())
            )/(rho1 + rho2),
            residualAlpha/2/L
        )
    );

    const volScalarField::Internal muI(rho1*nu1*rho2*nu2/(rho1*nu1 + rho2*nu2));

    const volScalarField::Internal limitedAlpha1
    (
        max(alpha1, interface_.phase1().residualAlpha())
    );

    const volScalarField::Internal limitedAlpha2
    (
        max(alpha2, interface_.phase2().residualAlpha())
    );

    const volScalarField::Internal muAlphaI
    (
        limitedAlpha1*rho1*nu1*limitedAlpha2*rho2*nu2
       /(limitedAlpha1*rho1*nu1 + limitedAlpha2*rho2*nu2)
    );

    const volScalarField::Internal ReI
    (
        (interface_.rho()()*interface_.magUr()()())
       /(magGradI*limitedAlpha1*limitedAlpha2*muI)
    );

    const volScalarField::Internal lambda(m_*ReI + n_*muAlphaI/muI);

    tmp<volScalarField> tK
    (
        volScalarField::New
        (
            "K",
            mesh,
            dimensionedScalar(dimK, 0),
            zeroGradientFvPatchField<scalar>::typeName
        )
    );

    tK.ref().ref() = lambda*sqr(magGradI)*muI;
    tK.ref().correctBoundaryConditions();

    return tK;
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::segregated::Kf() const
{
    return fvc::interpolate(K());
}


// ************************************************************************* //

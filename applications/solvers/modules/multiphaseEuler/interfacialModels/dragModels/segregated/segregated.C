/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

    const volScalarField& rho1(interface_.phase1().rho());
    const volScalarField& rho2(interface_.phase2().rho());

    tmp<volScalarField> tnu1(interface_.phase1().thermo().nu());
    tmp<volScalarField> tnu2(interface_.phase2().thermo().nu());

    const volScalarField& nu1(tnu1());
    const volScalarField& nu2(tnu2());

    volScalarField L
    (
        IOobject
        (
            "L",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimLength, 0),
        zeroGradientFvPatchField<scalar>::typeName
    );
    L.primitiveFieldRef() = cbrt(mesh.V());
    L.correctBoundaryConditions();

    const dimensionedScalar residualAlpha
    (
        (
            interface_.phase1().residualAlpha()
          + interface_.phase2().residualAlpha()
        )/2
    );

    const volScalarField I1(alpha1/max(alpha1 + alpha2, residualAlpha));
    const volScalarField I2(alpha2/max(alpha1 + alpha2, residualAlpha));
    const volScalarField magGradI
    (
        max
        (
            (rho2*mag(fvc::grad(I1)) + rho1*mag(fvc::grad(I2)))/(rho1 + rho2),
            residualAlpha/2/L
        )
    );

    const volScalarField muI(rho1*nu1*rho2*nu2/(rho1*nu1 + rho2*nu2));

    const volScalarField limitedAlpha1
    (
        max(alpha1, interface_.phase1().residualAlpha())
    );

    const volScalarField limitedAlpha2
    (
        max(alpha2, interface_.phase2().residualAlpha())
    );

    const volScalarField muAlphaI
    (
        limitedAlpha1*rho1*nu1*limitedAlpha2*rho2*nu2
       /(limitedAlpha1*rho1*nu1 + limitedAlpha2*rho2*nu2)
    );

    const volScalarField ReI
    (
        interface_.rho()*interface_.magUr()
       /(magGradI*limitedAlpha1*limitedAlpha2*muI)
    );

    const volScalarField lambda(m_*ReI + n_*muAlphaI/muI);

    return lambda*sqr(magGradI)*muI;
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::segregated::Kf() const
{
    return fvc::interpolate(K());
}


// ************************************************************************* //

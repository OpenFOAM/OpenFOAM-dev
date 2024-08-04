/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "JohnsonJacksonSchaefferFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(JohnsonJacksonSchaeffer, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        JohnsonJacksonSchaeffer,
        dictionary
    );
}
}
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::readCoeffs
(
    const dictionary& coeffDict
)
{
    Fr_.read(coeffDict);
    eta_.read(coeffDict);
    p_.read(coeffDict);

    phi_.read(coeffDict);
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::JohnsonJacksonSchaeffer
(
    const dictionary& coeffDict
)
:
    frictionalStressModel(coeffDict),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict),
    eta_("eta", dimless, coeffDict),
    p_("p", dimless, coeffDict),
    phi_("phi", dimless, coeffDict),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict)
{
    phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::~JohnsonJacksonSchaeffer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::frictionalPressure
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    const volScalarField& alpha = phase;

    return
        Fr_*pow(max(alpha - alphaMinFriction, scalar(0)), eta_)
       /pow(max(alphaMax - alpha, alphaDeltaMin_), p_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::frictionalPressurePrime
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    const volScalarField& alpha = phase;

    return Fr_*
    (
        eta_*pow(max(alpha - alphaMinFriction, scalar(0)), eta_ - 1)
       *(alphaMax - alpha)
      + p_*pow(max(alpha - alphaMinFriction, scalar(0)), eta_)
    )/pow(max(alphaMax - alpha, alphaDeltaMin_), p_ + 1);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::nu
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D
) const
{
    const volScalarField& alpha = phase;

    tmp<volScalarField> tnu
    (
        volScalarField::New
        (
            IOobject::groupName
            (
                Foam::typedName<frictionalStressModel>("nu"),
                phase.group()
            ),
            phase.mesh(),
            dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
        )
    );

    volScalarField& nuf = tnu.ref();

    forAll(D, celli)
    {
        if (alpha[celli] > alphaMinFriction.value())
        {
            nuf[celli] =
                0.5*pf[celli]*sin(phi_.value())
               /(
                    sqrt((1.0/3.0)*sqr(tr(D[celli])) - invariantII(D[celli]))
                  + small
                );
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            nufBf[patchi] =
                (
                    pf.boundaryField()[patchi]*sin(phi_.value())
                   /(
                        mag(phase.U()().boundaryField()[patchi].snGrad())
                      + small
                    )
                );
        }
    }

    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


// ************************************************************************* //

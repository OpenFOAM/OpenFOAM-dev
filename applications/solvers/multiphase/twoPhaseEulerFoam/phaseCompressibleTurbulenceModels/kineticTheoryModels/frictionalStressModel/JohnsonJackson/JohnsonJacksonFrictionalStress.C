/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "JohnsonJacksonFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(JohnsonJackson, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        JohnsonJackson,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
JohnsonJackson
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict_.lookup("Fr")),
    eta_("eta", dimless, coeffDict_.lookup("eta")),
    p_("p", dimless, coeffDict_.lookup("p")),
    phi_("phi", dimless, coeffDict_.lookup("phi"))
{
    phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
~JohnsonJackson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
frictionalPressure
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{

    return
        Fr_*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta_)
       /pow(max(alphaMax - alpha1, scalar(5.0e-2)), p_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
frictionalPressurePrime
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    return Fr_*
    (
        eta_*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta_ - 1.0)
       *(alphaMax-alpha1)
      + p_*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta_)
    )/pow(max(alphaMax - alpha1, scalar(5.0e-2)), p_ + 1.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::nu
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D
) const
{
    return dimensionedScalar("0.5", dimTime, 0.5)*pf*sin(phi_);
}


bool Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::read()
{
    coeffDict_ <<= dict_.subDict(typeName + "Coeffs");

    Fr_.read(coeffDict_);
    eta_.read(coeffDict_);
    p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;

    return true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "KochFriedlander.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{
namespace sinteringModels
{
    defineTypeNameAndDebug(KochFriedlander, 0);
    addToRunTimeSelectionTable(sinteringModel, KochFriedlander, dictionary);
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::sinteringModels::KochFriedlander::
KochFriedlander
(
    const dictionary& dict,
    const fractal& fractalShape
)
:
    sinteringModel(dict, fractalShape),
    dict_(dict.subDict(type() + "Coeffs")),
    Cs_(dict_.lookup<scalar>("Cs")),
    n_(dict_.lookup<scalar>("n")),
    m_(dict_.lookup<scalar>("m")),
    Ta_(dict_.lookup<scalar>("Ta"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::sinteringModels::KochFriedlander::
~KochFriedlander()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::shapeModels::sinteringModels::KochFriedlander::tau() const
{
    tmp<volScalarField::Internal> tTau
    (
        volScalarField::Internal::New
        (
            "tau",
            fractal_.SizeGroup().mesh(),
            dimensionedScalar(dimTime, 0.0)
        )
    );

    volScalarField::Internal& tau = tTau.ref();

    const sizeGroup& fi = fractal_.SizeGroup();
    const volScalarField& kappai = fractal_.fld();
    const volScalarField& T = fractal_.SizeGroup().phase().thermo().T();

    forAll(tau, celli)
    {
        tau[celli] =
            Cs_*pow(6.0/max(6.0/fi.dSph().value(), kappai[celli]), n_)
           *pow(T[celli], m_)*exp(Ta_/T[celli]);
    }

    return tTau;
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::shapeModels::sinteringModels::KochFriedlander::R() const
{
    const sizeGroup& fi = fractal_.SizeGroup();
    const volScalarField& kappai = fractal_.fld();
    const volScalarField& alpha = fi.phase();

    volScalarField::Internal R
    (
        IOobject
        (
            "KochFriedlander:R",
            fi.time().timeName(),
            fi.mesh()
        ),
        fi.mesh(),
        dimensionedScalar(inv(dimTime), 0)
    );

    volScalarField::Internal tau(this->tau());

    forAll(R, celli)
    {
        R[celli] = fi[celli]*alpha[celli]/tau[celli];
    }

    return fvm::Sp(R, kappai) - 6.0/fi.dSph()*R;
}


// ************************************************************************* //

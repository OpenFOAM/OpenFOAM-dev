/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SpalartAllmarasDDES.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmarasDDES, 0);
addToRunTimeSelectionTable(LESModel, SpalartAllmarasDDES, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasDDES::rd
(
    const volScalarField& visc,
    const volScalarField& S
) const
{
    return min
    (
        visc
        /(
            max
            (
                S,
                dimensionedScalar("SMALL", S.dimensions(), SMALL)
            )*sqr(kappa_*y_)
          + dimensionedScalar
            (
                "ROOTVSMALL",
                dimensionSet(0, 2 , -1, 0, 0),
                ROOTVSMALL
            )
        ),
        scalar(10)
    );
}


tmp<volScalarField> SpalartAllmarasDDES::fd(const volScalarField& S) const
{
    return 1 - tanh(pow3(8*rd(nuEff(), S)));
}


tmp<volScalarField> SpalartAllmarasDDES::S(const volTensorField& gradU) const
{
    return sqrt(2.0)*mag(symm(gradU));
}


tmp<volScalarField> SpalartAllmarasDDES::dTilda(const volScalarField& S) const
{
    return max
    (
        y_
      - fd(S)
       *max(y_ - CDES_*delta(), dimensionedScalar("zero", dimLength, 0)),
        dimensionedScalar("small", dimLength, SMALL)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmarasDDES::SpalartAllmarasDDES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    SpalartAllmaras(U, phi, transport, turbulenceModelName, modelName)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

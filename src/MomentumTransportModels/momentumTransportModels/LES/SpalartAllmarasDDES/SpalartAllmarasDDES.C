/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasDDES<BasicMomentumTransportModel>::rd
(
    const volInternalScalarField& magGradU
) const
{
    return volInternalScalarField::New
    (
        typedName("rd"),
        min
        (
            this->nuEff()()
           /(
                max
                (
                    magGradU,
                    dimensionedScalar(magGradU.dimensions(), small)
                )
               *sqr(this->kappa_*this->y()())
            ),
            scalar(10)
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasDDES<BasicMomentumTransportModel>::fd
(
    const volInternalScalarField& magGradU
) const
{
    return volInternalScalarField::New
    (
        typedName("fd"),
        1 - tanh(pow3(8*rd(magGradU)))
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasDDES<BasicMomentumTransportModel>::dTilda
(
    const volInternalScalarField& chi,
    const volInternalScalarField& fv1,
    const volInternalTensorField& gradU
) const
{
    return volInternalScalarField::New
    (
        typedName("dTilda"),
        max
        (
            this->y()
          - fd(mag(gradU))
           *max
            (
                this->y()() - this->CDES_*this->delta()(),
                dimensionedScalar(dimensions::length, 0)
            ),
            dimensionedScalar(dimensions::length, small)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
SpalartAllmarasDDES<BasicMomentumTransportModel>::SpalartAllmarasDDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    SpalartAllmarasDES<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //

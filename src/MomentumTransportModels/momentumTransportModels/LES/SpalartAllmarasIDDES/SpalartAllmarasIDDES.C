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

#include "SpalartAllmarasIDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::IDDESalpha() const
{
    return volInternalScalarField::New
    (
        typedName("alpha"),
        max(0.25 - this->y()()/IDDESDelta_.hmax(), scalar(-5))
    );
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::ft
(
    const volInternalScalarField& magGradU
) const
{
    return volInternalScalarField::New
    (
        typedName("ft"),
        tanh(pow3(sqr(ct_)*rd(this->nut_, magGradU)))
    );
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::fl
(
    const volInternalScalarField& magGradU
) const
{
    return volInternalScalarField::New
    (
        typedName("fl"),
        tanh(pow(sqr(cl_)*rd(this->nu(), magGradU), 10))
    );
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::rd
(
    const volInternalScalarField& nur,
    const volInternalScalarField& magGradU
) const
{
    return volInternalScalarField::New
    (
        typedName("rd"),
        min
        (
            nur
           /(
               max
               (
                   magGradU,
                   dimensionedScalar(magGradU.dimensions(), small)
               )*sqr(this->kappa_*this->y()())
            ),
            scalar(10)
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::fd
(
    const volInternalScalarField& magGradU
) const
{
    return volInternalScalarField::New
    (
        typedName("fd"),
        1 - tanh(pow3(8*rd(this->nuEff(), magGradU)))
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::dTilda
(
    const volInternalScalarField& chi,
    const volInternalScalarField& fv1,
    const volInternalTensorField& gradU
) const
{
    const volInternalScalarField alpha(IDDESalpha());

    const volInternalScalarField expTerm
    (
        typedName("expTerm"),
        exp(sqr(alpha))
    );

    const volInternalScalarField magGradU(typedName("magGradU"), mag(gradU));

    tmp<volInternalScalarField> fHill
    (
        volInternalScalarField::New
        (
            typedName("fHill"),
            2*(pos0(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0))
        )
    );

    tmp<volInternalScalarField> fStep
    (
        volInternalScalarField::New
        (
            typedName("fStep"),
            min(2*pow(expTerm, -9.0), scalar(1))
        )
    );

    const volInternalScalarField fHyb
    (
        typedName("fHyb"),
        max(1 - fd(magGradU), fStep)
    );

    tmp<volInternalScalarField> fAmp
    (
        volInternalScalarField::New
        (
            typedName("fAmp"),
            1 - max(ft(magGradU), fl(magGradU))
        )
    );

    tmp<volInternalScalarField> fRestore
    (
        volInternalScalarField::New
        (
            typedName("fRestore"),
            max(fHill - 1, scalar(0))*fAmp
        )
    );

    // IGNORING ft2 terms
    const volInternalScalarField Psi
    (
        typedName("Psi"),
        sqrt
        (
            min
            (
                scalar(100),
                (
                    1
                  - this->Cb1_*this->fv2(chi, fv1)
                   /(this->Cw1_*sqr(this->kappa_)*fwStar_)
                 )/max(small, fv1)
            )
        )
    );

    return volInternalScalarField::New
    (
        typedName("dTilda"),
        max
        (
            dimensionedScalar(dimensions::length, small),
            fHyb*(1 + fRestore*Psi)*this->y()
          + (1 - fHyb)*this->CDES_*Psi*this->delta()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::SpalartAllmarasIDDES
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
    ),
    fwStar_("fwStar", this->typeDict(type), 0.424),
    cl_("cl", this->typeDict(type), 3.55),
    ct_("ct", this->typeDict(type), 1.63),
    IDDESDelta_(refCast<IDDESDelta>(this->delta_()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool SpalartAllmarasIDDES<BasicMomentumTransportModel>::read()
{
    if (SpalartAllmarasDES<BasicMomentumTransportModel>::read())
    {
        fwStar_.readIfPresent(this->typeDict());
        cl_.readIfPresent(this->typeDict());
        ct_.readIfPresent(this->typeDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //

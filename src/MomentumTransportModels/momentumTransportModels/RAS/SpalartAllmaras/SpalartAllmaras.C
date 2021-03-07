/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "SpalartAllmaras.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> SpalartAllmaras<BasicMomentumTransportModel>::chi() const
{
    return volScalarField::New(modelName("chi"), nuTilda_/this->nu());
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> SpalartAllmaras<BasicMomentumTransportModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(modelName("chi3"), pow3(chi));
    return volScalarField::New(modelName("fv1"), chi3/(chi3 + pow3(Cv1_)));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> SpalartAllmaras<BasicMomentumTransportModel>::fv2
(
    const volScalarField::Internal& chi,
    const volScalarField::Internal& fv1
) const
{
    return volScalarField::Internal::New
    (
        modelName("fv2"),
        1.0 - chi/(1.0 + chi*fv1)
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmaras<BasicMomentumTransportModel>::Stilda
(
    const volScalarField::Internal& chi,
    const volScalarField::Internal& fv1
) const
{
    const volScalarField::Internal Omega
    (
        modelName("Omega"),
        ::sqrt(2.0)*mag(skew(fvc::grad(this->U_)().v()))
    );

    return volScalarField::Internal::New
    (
        modelName("Stilda"),
        (
            max
            (
                Omega
              + fv2(chi, fv1)*nuTilda_/sqr(kappa_*y_),
                Cs_*Omega
            )
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> SpalartAllmaras<BasicMomentumTransportModel>::fw
(
    const volScalarField::Internal& Stilda
) const
{
    const volScalarField::Internal r
    (
        modelName("r"),
        min
        (
            nuTilda_()
           /(
               max
               (
                   Stilda,
                   dimensionedScalar(Stilda.dimensions(), small)
               )
              *sqr(kappa_*y_)
            ),
            scalar(10.0)
        )
    );

    const volScalarField::Internal g(modelName("g"), r + Cw2_*(pow6(r) - r));

    return volScalarField::Internal::New
    (
        modelName("fw"),
        g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0)
    );
}


template<class BasicMomentumTransportModel>
void SpalartAllmaras<BasicMomentumTransportModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = nuTilda_*fv1;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void SpalartAllmaras<BasicMomentumTransportModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
SpalartAllmaras<BasicMomentumTransportModel>::SpalartAllmaras
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool SpalartAllmaras<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
SpalartAllmaras<BasicMomentumTransportModel>::DnuTildaEff() const
{
    return volScalarField::New
    (
        "DnuTildaEff",
        (nuTilda_ + this->nu())/sigmaNut_
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> SpalartAllmaras<BasicMomentumTransportModel>::k() const
{
    return volScalarField::New
    (
        "k",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -2, 0, 0), 0)
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
SpalartAllmaras<BasicMomentumTransportModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return volScalarField::New
    (
        "epsilon",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -3, 0, 0), 0)
    );
}


template<class BasicMomentumTransportModel>
void SpalartAllmaras<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    const volScalarField::Internal Stilda(this->Stilda(chi, fv1));

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1_*alpha()*rho()*Stilda*nuTilda_()
      - fvm::Sp(Cw1_*alpha()*rho()*fw(Stilda)*nuTilda_()/sqr(y_), nuTilda_)
      + fvModels.source(alpha, rho, nuTilda_)
    );

    nuTildaEqn.ref().relax();
    fvConstraints.constrain(nuTildaEqn.ref());
    solve(nuTildaEqn);
    fvConstraints.constrain(nuTilda_);
    bound(nuTilda_, dimensionedScalar(nuTilda_.dimensions(), 0));
    nuTilda_.correctBoundaryConditions();

    correctNut(fv1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

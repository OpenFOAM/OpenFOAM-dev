/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(LESModel, SpalartAllmaras, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void SpalartAllmaras::updateSubGridScaleFields()
{
    nuSgs_.internalField() = fv1()*nuTilda_.internalField();
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::fv1() const
{
    const volScalarField chi3("chi3", pow3(nuTilda_/nu()));
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SpalartAllmaras::fv2() const
{
    if (ashfordCorrection_)
    {
        return 1/pow3(scalar(1) + nuTilda_/(Cv2_*nu()));
    }
    else
    {
        const volScalarField chi("chi", nuTilda_/nu());
        return 1.0 - chi/(1.0 + chi*fv1());
    }
}


tmp<volScalarField> SpalartAllmaras::fv3() const
{
    if (ashfordCorrection_)
    {
        const volScalarField chi("chi", nuTilda_/nu());
        const volScalarField chiByCv2(chi/Cv2_);

        return
            (scalar(1) + chi*fv1())
           *(1/Cv2_)
           *(3*(scalar(1) + chiByCv2) + sqr(chiByCv2))
           /pow3(scalar(1) + chiByCv2);
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fv3",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("fv3", dimless, 1),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}


tmp<volScalarField> SpalartAllmaras::S(const volTensorField& gradU) const
{
    return sqrt(2.0)*mag(skew(gradU));
}


tmp<volScalarField> SpalartAllmaras::STilda
(
    const volScalarField& S,
    const volScalarField& dTilda
) const
{
    return fv3()*S + fv2()*nuTilda_/sqr(kappa_*dTilda);
}


tmp<volScalarField> SpalartAllmaras::r
(
    const volScalarField& visc,
    const volScalarField& S,
    const volScalarField& dTilda
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
           )
          *sqr(kappa_*dTilda)
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


tmp<volScalarField> SpalartAllmaras::fw
(
    const volScalarField& S,
    const volScalarField& dTilda
) const
{
    const volScalarField r(this->r(nuTilda_, S, dTilda));
    const volScalarField g(r + Cw2_*(pow6(r) - r));

    return g*pow((1 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


tmp<volScalarField> SpalartAllmaras::dTilda(const volScalarField&) const
{
    return min(CDES_*delta(), y_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmaras::SpalartAllmaras
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    DESModel(modelName, U, phi, transport, turbulenceModelName),

    sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            coeffDict_,
            0.622
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            coeffDict_,
            7.1
        )
    ),
    Cv2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv2",
            coeffDict_,
            5.0
        )
    ),
    CDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES",
            coeffDict_,
            0.65
        )
    ),
    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.07
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            2.0
        )
    ),

    ashfordCorrection_(coeffDict_.lookupOrDefault("ashfordCorrection", true)),

    y_(mesh_),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    updateSubGridScaleFields();

    if (ashfordCorrection_)
    {
        Info<< "    Employing Ashford correction" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SpalartAllmaras::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    if (mesh_.changing())
    {
        y_.correct();
        y_.boundaryField() = max(y_.boundaryField(), VSMALL);
    }

    const volScalarField S(this->S(gradU));
    const volScalarField dTilda(this->dTilda(S));
    const volScalarField STilda(this->STilda(S, dTilda));

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(nuTilda_)
      + fvm::div(phi(), nuTilda_)
      - fvm::laplacian
        (
            (nuTilda_ + nu())/sigmaNut_,
            nuTilda_,
            "laplacian(DnuTildaEff,nuTilda)"
        )
      - Cb2_/sigmaNut_*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1_*STilda*nuTilda_
      - fvm::Sp(Cw1_*fw(STilda, dTilda)*nuTilda_/sqr(dTilda), nuTilda_)
    );

    nuTildaEqn().relax();
    nuTildaEqn().solve();

    bound(nuTilda_, dimensionedScalar("zero", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    updateSubGridScaleFields();
}


tmp<volScalarField> SpalartAllmaras::k() const
{
    return sqr(nuSgs()/ck_/dTilda(S(fvc::grad(U()))));
}


tmp<volScalarField> SpalartAllmaras::epsilon() const
{
    return 2*nuEff()*magSqr(symm(fvc::grad(U())));
}


tmp<volSymmTensorField> SpalartAllmaras::B() const
{
    return ((2.0/3.0)*I)*k() - nuSgs()*twoSymm(fvc::grad(U()));
}


tmp<volSymmTensorField> SpalartAllmaras::devReff() const
{
    return -nuEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<fvVectorMatrix> SpalartAllmaras::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> SpalartAllmaras::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool SpalartAllmaras::read()
{
    if (LESModel::read())
    {
        sigmaNut_.readIfPresent(coeffDict());
        kappa_.readIfPresent(*this);

        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        Cv2_.readIfPresent(coeffDict());
        CDES_.readIfPresent(coeffDict());
        ck_.readIfPresent(coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());

        ashfordCorrection_.readIfPresent("ashfordCorrection", coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


tmp<volScalarField> SpalartAllmaras::LESRegion() const
{
    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neg(dTilda(S(fvc::grad(U_))) - y_)
//            mesh_,
//            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    return tLESRegion;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

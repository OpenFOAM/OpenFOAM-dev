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
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(LESModel, SpalartAllmaras, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void SpalartAllmaras::updateSubGridScaleFields()
{
    muSgs_.internalField() = rho()*fv1()*nuTilda_.internalField();
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


tmp<volScalarField> SpalartAllmaras::fv1() const
{
    volScalarField chi3(pow3(rho()*nuTilda_/mu()));
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SpalartAllmaras::fv2() const
{
    if (ashfordCorrection_)
    {
        return 1.0/pow3(scalar(1) + rho()*nuTilda_/(mu()*Cv2_));
    }
    else
    {
        const volScalarField chi("chi", rho()*nuTilda_/mu());
        return 1.0 - chi/(1.0 + chi*fv1());
    }
}


tmp<volScalarField> SpalartAllmaras::fv3() const
{
    if (ashfordCorrection_)
    {
        volScalarField chi(rho()*nuTilda_/mu());
        volScalarField chiByCv2((1/Cv2_)*chi);

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


tmp<volScalarField> SpalartAllmaras::fw(const volScalarField& Stilda) const
{
    volScalarField r
    (
        min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa_*dTilda_)
            ),
            scalar(10.0)
        )
    );
    r.boundaryField() == 0.0;

    volScalarField g(r + Cw2_*(pow6(r) - r));

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmaras::SpalartAllmaras
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermoPhysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    DESModel(modelName, rho, U, phi, thermoPhysicalModel, turbulenceModelName),

    sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
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
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            *this,
            0.41
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

    dTilda_(min(CDES_*delta(), wallDist(mesh_).y())),
    muSgs_
    (
        IOobject
        (
            "muSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    alphaSgs_
    (
        IOobject
        (
            "alphaSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    updateSubGridScaleFields();

    printCoeffs();

    if (ashfordCorrection_)
    {
        Info<< "    Employing Ashford correction" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> SpalartAllmaras::B() const
{
    return ((2.0/3.0)*I)*k() - (muSgs_/rho())*dev(twoSymm(fvc::grad(U())));
}


tmp<volSymmTensorField> SpalartAllmaras::devRhoBeff() const
{
    return -muEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<volScalarField> SpalartAllmaras::epsilon() const
{
    return 2*muEff()/rho()*magSqr(symm(fvc::grad(U())));
}


tmp<fvVectorMatrix> SpalartAllmaras::divDevRhoBeff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U) - fvc::div(muEff()*dev2(T(fvc::grad(U))))
    );
}


void SpalartAllmaras::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();
    LESModel::correct(gradU);

    if (mesh_.changing())
    {
        dTilda_ = min(CDES_*delta(), wallDist(mesh_).y());
    }

    volScalarField Stilda
    (
        fv3()*::sqrt(2.0)*mag(skew(gradU)) + fv2()*nuTilda_/sqr(kappa_*dTilda_)
    );

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(rho(), nuTilda_)
      + fvm::div(phi(), nuTilda_)
      - fvm::laplacian
        (
            (nuTilda_*rho() + mu())/sigmaNut_,
            nuTilda_,
            "laplacian(DnuTildaEff,nuTilda)"
        )
      - rho()*Cb2_/sigmaNut_*magSqr(fvc::grad(nuTilda_))
     ==
        rho()*Cb1_*Stilda*nuTilda_
      - fvm::Sp(rho()*Cw1_*fw(Stilda)*nuTilda_/sqr(dTilda_), nuTilda_)
    );

    nuTildaEqn().relax();
    nuTildaEqn().solve();

    bound(nuTilda_, dimensionedScalar("zero", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    updateSubGridScaleFields();
}


bool SpalartAllmaras::read()
{
    if (LESModel::read())
    {
        sigmaNut_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());

        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        Cv2_.readIfPresent(coeffDict());
        CDES_.readIfPresent(coeffDict());
        ck_.readIfPresent(coeffDict());
        kappa_.readIfPresent(*this);
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
    volScalarField wd(wallDist(mesh_).y());

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
            neg(min(CDES_*delta(), wd) - wd)
//            mesh_,
//            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    return tLESRegion;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "DarcyForchheimer.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(DarcyForchheimer, 0);
        addToRunTimeSelectionTable(porosityModel, DarcyForchheimer, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::DarcyForchheimer
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& coeffDict,
    const word& cellZoneName
)
:
    porosityModel(name, mesh, dict, coeffDict, cellZoneName),
    dXYZ_("d", dimless/sqr(dimLength), coeffDict),
    fXYZ_("f", dimless/dimLength, coeffDict),
    rhoName_(coeffDict.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffDict.lookupOrDefault<word>("mu", "mu")),
    nuName_(coeffDict.lookupOrDefault<word>("nu", "nu"))
{
    adjustNegativeResistance(dXYZ_);
    adjustNegativeResistance(fXYZ_);

    calcTransformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::~DarcyForchheimer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::DarcyForchheimer::calcTransformModelData()
{
    if (coordSys_.R().uniform())
    {
        D_.setSize(1);
        F_.setSize(1);

        D_[0] = Zero;
        D_[0].xx() = dXYZ_.value().x();
        D_[0].yy() = dXYZ_.value().y();
        D_[0].zz() = dXYZ_.value().z();

        D_[0] = coordSys_.R().transform(Zero, D_[0]);

        // leading 0.5 is from 1/2*rho
        F_[0] = Zero;
        F_[0].xx() = 0.5*fXYZ_.value().x();
        F_[0].yy() = 0.5*fXYZ_.value().y();
        F_[0].zz() = 0.5*fXYZ_.value().z();

        F_[0] = coordSys_.R().transform(Zero, F_[0]);
    }
    else
    {
        const labelList& cells = mesh_.cellZones()[zoneName_];

        D_.setSize(cells.size());
        F_.setSize(cells.size());

        forAll(cells, i)
        {
            D_[i] = Zero;
            D_[i].xx() = dXYZ_.value().x();
            D_[i].yy() = dXYZ_.value().y();
            D_[i].zz() = dXYZ_.value().z();

            // leading 0.5 is from 1/2*rho
            F_[i] = Zero;
            F_[i].xx() = 0.5*fXYZ_.value().x();
            F_[i].yy() = 0.5*fXYZ_.value().y();
            F_[i].zz() = 0.5*fXYZ_.value().z();
        }

        const coordinateRotation& R = coordSys_.R
        (
            UIndirectList<vector>(mesh_.C(), cells)()
        );

        D_ = R.transform(D_);
        F_ = R.transform(F_);
    }

    if (debug && (mesh_.time().writeTime() || mesh_.time().timeIndex() == 0))
    {
        volTensorField Dout
        (
            IOobject
            (
                typedName("D"),
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(dXYZ_.dimensions(), Zero)
        );
        volTensorField Fout
        (
            IOobject
            (
                typedName("F"),
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(fXYZ_.dimensions(), Zero)
        );

        UIndirectList<tensor>(Dout, mesh_.cellZones()[zoneName_]) = D_;
        UIndirectList<tensor>(Fout, mesh_.cellZones()[zoneName_]) = F_;

        Dout.write();
        Fout.write();
    }
}


void Foam::porosityModels::DarcyForchheimer::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu, U);
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu, U);
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho, U);
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);
        const volScalarField& mu = mesh_.lookupObject<volScalarField>(muName);

        apply(AU, rho, mu, U);
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(AU, geometricOneField(), nu, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(AU, geometricOneField(), mu/rho, U);
        }
    }
}


// ************************************************************************* //

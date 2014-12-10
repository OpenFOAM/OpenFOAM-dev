/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "fvMatrices.H"
#include "volFields.H"
#include "DimensionedField.H"

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh> >
Foam::coalCloudList::UTrans() const
{
    tmp<DimensionedField<vector, volMesh> > tfld
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                "UTransEff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimMass*dimVelocity, vector::zero)
        )
    );

    DimensionedField<vector, volMesh>& fld = tfld();

    forAll(*this, i)
    {
        fld += operator[](i).UTrans();
    }

    return tfld;
}


Foam::tmp<Foam::fvVectorMatrix> Foam::coalCloudList::SU
(
    volVectorField& U
) const
{
    tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
    fvVectorMatrix& fvm = tfvm();

    forAll(*this, i)
    {
        fvm += operator[](i).SU(U);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::coalCloudList::hsTrans() const
{
    tmp<DimensionedField<scalar, volMesh> > tfld
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "hsTransEff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& fld = tfld();

    forAll(*this, i)
    {
        fld += operator[](i).hsTrans();
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::Sh
(
    volScalarField& hs
) const
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
    fvScalarMatrix& fvm = tfvm();

    forAll(*this, i)
    {
        fvm += operator[](i).Sh(hs);
    }

    return tfvm;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::SYi
(
    const label ii,
    volScalarField& Yi
) const
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass/dimTime));
    fvScalarMatrix& fvm = tfvm();

    forAll(*this, i)
    {
        fvm += operator[](i).SYi(ii, Yi);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::coalCloudList::rhoTrans() const
{
    tmp<DimensionedField<scalar, volMesh> > tfld
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "rhoTransEff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& fld = tfld();

    forAll(*this, i)
    {
        forAll(operator[](i).rhoTrans(), j)
        {
            fld += operator[](i).rhoTrans()[j];
        }
    }

    return tfld;
}




Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::coalCloudList::Srho() const
{
    tmp<DimensionedField<scalar, volMesh> > tfld
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "rhoTransEff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& fld = tfld();

    forAll(*this, i)
    {
        fld += operator[](i).Srho();
    }

    return tfld;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::coalCloudList::Srho
(
    const label i
) const
{
    tmp<DimensionedField<scalar, volMesh> > tfld
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "rhoTransEff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& fld = tfld();

    forAll(*this, j)
    {
        fld += operator[](j).Srho(i);
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::Srho
(
    volScalarField& rho
) const
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
    fvScalarMatrix& fvm = tfvm();

    forAll(*this, i)
    {
        fvm += operator[](i).Srho(rho);
    }

    return tfvm;
}


// ************************************************************************* //

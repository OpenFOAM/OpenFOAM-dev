/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "parcelCloudList.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "wordIOList.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::parcelCloudList::defaultCloudName("cloud");

const Foam::wordList Foam::parcelCloudList::defaultCloudNames(1, "cloud");

const Foam::word Foam::parcelCloudList::cloudNamesName("clouds");


// * * * * * * * * * * Private Static Member Functions * * * * * * * * * * * //

Foam::wordList Foam::parcelCloudList::cloudNames(const objectRegistry& db)
{
    typeIOobject<wordGlobalIOList> cloudNamesIO
    (
        cloudNamesName,
        db.time().constant(),
        db,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (cloudNamesIO.headerOk())
    {
        return wordGlobalIOList(cloudNamesIO);
    }

    typeIOobject<IOdictionary> cloudIO
    (
        defaultCloudName + "Properties",
        db.time().constant(),
        db,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (cloudIO.headerOk())
    {
        return wordList(1, defaultCloudName);
    }

    FatalErrorInFunction
        << "Cloud properties were not found in either "
        << cloudNamesIO.relativeObjectPath() << " or "
        << cloudIO.relativeObjectPath() << exit(FatalError);

    return wordList::null();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parcelCloudList::parcelCloudList
(
    const wordList& cloudNames,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g
)
:
    PtrList<parcelCloud>(),
    mesh_(rho.mesh())
{
    this->setSize(cloudNames.size());

    forAll(cloudNames, i)
    {
        this->set(i, parcelCloud::New(cloudNames[i], rho, U, mu, g));
    }
}


Foam::parcelCloudList::parcelCloudList
(
    const wordList& cloudNames,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo
)
:
    PtrList<parcelCloud>(),
    mesh_(rho.mesh())
{
    this->setSize(cloudNames.size());

    forAll(cloudNames, i)
    {
        this->set(i, parcelCloud::New(cloudNames[i], rho, U, g, carrierThermo));
    }
}


Foam::parcelCloudList::parcelCloudList
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g
)
:
    parcelCloudList(cloudNames(rho.mesh()), rho, U, mu, g)
{}


Foam::parcelCloudList::parcelCloudList
(
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo
)
:
    parcelCloudList(cloudNames(rho.mesh()), rho, U, g, carrierThermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parcelCloudList::~parcelCloudList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::tmp<Foam::volScalarField> Foam::parcelCloudList::theta() const
{
    tmp<volScalarField> ttheta
    (
        volScalarField::New
        (
            cloudNamesName + ":theta",
            mesh_,
            dimensionedScalar(dimless, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    forAll(*this, i)
    {
        ttheta.ref() += operator[](i).theta();
    }
    return ttheta;
}


Foam::tmp<Foam::fvVectorMatrix> Foam::parcelCloudList::SU
(
    const volVectorField& U
) const
{
    tmp<fvVectorMatrix> tSU(new fvVectorMatrix(U, dimMass*dimAcceleration));
    forAll(*this, i)
    {
        tSU.ref() += operator[](i).SU(U);
    }
    return tSU;
}


Foam::tmp<Foam::volVectorField::Internal> Foam::parcelCloudList::UTrans() const
{
    tmp<volVectorField::Internal> tUTrans
    (
        volVectorField::Internal::New
        (
            cloudNamesName + ":UTrans",
            mesh_,
            dimensionedVector(dimMass*dimVelocity, Zero)
        )
    );
    forAll(*this, i)
    {
        tUTrans.ref() += operator[](i).UTrans();
    }
    return tUTrans;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::parcelCloudList::UCoeff() const
{
    tmp<volScalarField::Internal> tUCoeff
    (
        volScalarField::Internal::New
        (
            cloudNamesName + ":UCoeff",
            mesh_,
            dimensionedScalar(dimMass, Zero)
        )
    );
    forAll(*this, i)
    {
        tUCoeff.ref() += operator[](i).UCoeff();
    }
    return tUCoeff;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::parcelCloudList::Sh
(
    const volScalarField& hs
) const
{
    tmp<fvScalarMatrix> tSh(new fvScalarMatrix(hs, dimEnergy/dimTime));
    forAll(*this, i)
    {
        tSh.ref() += operator[](i).Sh(hs);
    }
    return tSh;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::parcelCloudList::hsTrans() const
{
    tmp<volScalarField::Internal> thsTrans
    (
        volScalarField::Internal::New
        (
            cloudNamesName + ":hsTrans",
            mesh_,
            dimensionedScalar(dimEnergy, Zero)
        )
    );
    forAll(*this, i)
    {
        thsTrans.ref() += operator[](i).hsTrans();
    }
    return thsTrans;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::parcelCloudList::hsCoeff() const
{
    tmp<volScalarField::Internal> thsCoeff
    (
        volScalarField::Internal::New
        (
            cloudNamesName + ":hsCoeff",
            mesh_,
            dimensionedScalar(dimEnergy/dimTemperature, Zero)
        )
    );
    forAll(*this, i)
    {
        thsCoeff.ref() += operator[](i).hsCoeff();
    }
    return thsCoeff;
}


Foam::tmp<Foam::volScalarField> Foam::parcelCloudList::Ep() const
{
    tmp<volScalarField> tEp
    (
        volScalarField::New
        (
            cloudNamesName + ":radiation:Ep",
            mesh_,
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
        )
    );
    forAll(*this, i)
    {
        tEp.ref() += operator[](i).Ep();
    }
    return tEp;
}


Foam::tmp<Foam::volScalarField> Foam::parcelCloudList::ap() const
{
    tmp<volScalarField> tap
    (
        volScalarField::New
        (
            cloudNamesName + ":radiation:ap",
            mesh_,
            dimensionedScalar(dimless/dimLength, 0)
        )
    );
    forAll(*this, i)
    {
        tap.ref() += operator[](i).ap();
    }
    return tap;

}


Foam::tmp<Foam::volScalarField> Foam::parcelCloudList::sigmap() const
{
    tmp<volScalarField> tsigmap
    (
        volScalarField::New
        (
            cloudNamesName + ":radiation:sigmap",
            mesh_,
            dimensionedScalar(dimless/dimLength, 0)
        )
    );
    forAll(*this, i)
    {
        tsigmap.ref() += operator[](i).sigmap();
    }
    return tsigmap;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::parcelCloudList::SYi
(
    const label speciei,
    const volScalarField& Yi
) const
{
    tmp<fvScalarMatrix> tSYi(new fvScalarMatrix(Yi, dimMass/dimTime));
    forAll(*this, i)
    {
        tSYi.ref() += operator[](i).SYi(speciei, Yi);
    }
    return tSYi;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::parcelCloudList::Srho
(
    const volScalarField& rho
) const
{
    tmp<fvScalarMatrix> tSrho(new fvScalarMatrix(rho, dimMass/dimTime));
    forAll(*this, i)
    {
        tSrho.ref() += operator[](i).Srho(rho);
    }
    return tSrho;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::parcelCloudList::Srho() const
{
    tmp<volScalarField::Internal> tSrho
    (
        volScalarField::Internal::New
        (
            cloudNamesName + ":Srho",
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );
    forAll(*this, i)
    {
        tSrho.ref() += operator[](i).Srho();
    }
    return tSrho;
}


void Foam::parcelCloudList::info()
{
    forAll(*this, i)
    {
        operator[](i).info();
    }
}


void Foam::parcelCloudList::evolve()
{
    forAll(*this, i)
    {
        operator[](i).evolve();
    }
}


void Foam::parcelCloudList::storeGlobalPositions()
{
    forAll(*this, i)
    {
        operator[](i).storeGlobalPositions();
    }
}


void Foam::parcelCloudList::topoChange(const polyTopoChangeMap& map)
{
    forAll(*this, i)
    {
        operator[](i).topoChange(map);
    }
}


void Foam::parcelCloudList::mapMesh(const polyMeshMap& map)
{
    forAll(*this, i)
    {
        operator[](i).mapMesh(map);
    }
}


void Foam::parcelCloudList::distribute(const polyDistributionMap& map)
{
    forAll(*this, i)
    {
        operator[](i).distribute(map);
    }
}


// ************************************************************************* //

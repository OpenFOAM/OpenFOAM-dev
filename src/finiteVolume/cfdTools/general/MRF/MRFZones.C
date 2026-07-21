/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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

#include "MRFZones.H"
#include "volFields.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MRFZones, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::MRFZones::createIOobject(const fvMesh& mesh)
{
    typeIOobject<IOdictionary> io
    (
        "MRFProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        Info<< indentOrNl
            << "Constructing MRF zones from " << io.name()
            << endl;

        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZones::MRFZones(const fvMesh& mesh)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        TopoChangeableMeshObject,
        MRFZones,
        IOdictionary
    >
    (
        createIOobject(mesh),
        mesh
    ),
    PtrList<MRFZone>()
{
    reset(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFZones::reset(const dictionary& dict)
{
    printDictionary print(dict);

    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    if (count)
    {
        Info<< indentOrNl << "MRF zone list" << endl;
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            printDictionary print(modelDict);
            Info<< indent << "Constructing MRF zone " << name << endl;
            PtrList<MRFZone>::set
            (
                i++,
                new MRFZone(name, mesh(), modelDict)
            );
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::MRFZones::DDt
(
    const volVectorField& U
) const
{
    tmp<volVectorField> tDDt
    (
        volVectorField::New
        (
            "MRFZones:DDt",
            U.mesh(),
            dimensionedVector(U.dimensions()/dimensions::time, Zero)
        )
    );
    volVectorField& DDt = tDDt.ref();

    forAll(*this, i)
    {
        operator[](i).addCoriolis(U, DDt);
    }

    return tDDt;
}


Foam::tmp<Foam::volVectorField> Foam::MRFZones::DDt
(
    const volScalarField& rho,
    const volVectorField& U
) const
{
    return rho*DDt(U);
}


Foam::tmp<Foam::volVectorField>
Foam::MRFZones::centrifugalAcceleration() const
{
    tmp<volVectorField> tcentrifugalAcceleration
    (
        volVectorField::New
        (
            "MRFZones:centrifugalAcceleration",
            mesh(),
            dimensionedVector(dimensions::acceleration, Zero)
        )
    );
    volVectorField& centrifugalAcceleration = tcentrifugalAcceleration.ref();

    forAll(*this, i)
    {
        operator[](i).addCentrifugalAcceleration(centrifugalAcceleration);
    }

    return tcentrifugalAcceleration;
}



void Foam::MRFZones::makeRelative(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(U);
    }
}


void Foam::MRFZones::makeRelative(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZones::relative
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            Foam::New
            (
                tphi,
                "relative(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeRelative(rphi.ref());

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::FieldField<Foam::surfaceMesh::PatchField, Foam::scalar>>
Foam::MRFZones::relative
(
    const tmp<FieldField<surfaceMesh::PatchField, scalar>>& tphi
) const
{
    if (size())
    {
        tmp<FieldField<surfaceMesh::PatchField, scalar>> rphi
        (
            Foam::New(tphi, true)
        );

        forAll(*this, i)
        {
            operator[](i).makeRelative(rphi.ref());
        }

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<FieldField<surfaceMesh::PatchField, scalar>>(tphi, true);
    }
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::MRFZones::relative
(
    const tmp<Field<scalar>>& tphi,
    const label patchi
) const
{
    if (size())
    {
        tmp<Field<scalar>> rphi(Foam::New(tphi, true));

        forAll(*this, i)
        {
            operator[](i).makeRelative(rphi.ref(), patchi);
        }

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<Field<scalar>>(tphi, true);
    }
}


void Foam::MRFZones::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(rho, phi);
    }
}


void Foam::MRFZones::makeAbsolute(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(U);
    }
}


void Foam::MRFZones::makeAbsolute(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZones::absolute
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            Foam::New
            (
                tphi,
                "absolute(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeAbsolute(rphi.ref());

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


void Foam::MRFZones::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(rho, phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZones::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            Foam::New
            (
                tphi,
                "absolute(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeAbsolute(fvc::interpolate(rho), rphi.ref());

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


bool Foam::MRFZones::movePoints()
{
    return true;
}


void Foam::MRFZones::topoChange(const polyTopoChangeMap& map)
{
    forAll(*this, i)
    {
        operator[](i).topoChange(map);
    }
}


void Foam::MRFZones::mapMesh(const polyMeshMap& map)
{
    forAll(*this, i)
    {
        operator[](i).mapMesh(map);
    }
}


void Foam::MRFZones::distribute(const polyDistributionMap& map)
{
    forAll(*this, i)
    {
        operator[](i).distribute(map);
    }
}


bool Foam::MRFZones::read()
{
    if (regIOobject::read())
    {
        bool allOk = true;
        forAll(*this, i)
        {
            MRFZone& pm = this->operator[](i);
            bool ok = pm.read(this->subDict(pm.name()));
            allOk = (allOk && ok);
        }
        return allOk;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

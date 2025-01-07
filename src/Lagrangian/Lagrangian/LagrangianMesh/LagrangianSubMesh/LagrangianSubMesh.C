/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianSubMesh.H"
#include "LagrangianSubFields.H"
#include "LagrangianMesh.H"
#include "tracking.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LagrangianSubMesh, 0);
}


// * * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * //

Foam::LagrangianSubMesh::LagrangianSubMesh
(
    const LagrangianMesh& mesh,
    const LagrangianGroup group,
    const label size,
    const label start,
    const label index
)
:
    GeoMesh<polyMesh>(mesh.mesh()),
    mesh_(mesh),
    group_(group),
    size_(size),
    start_(start),
    index_(index)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianSubMesh::LagrangianSubMesh
(
    const LagrangianMesh& mesh,
    const LagrangianGroup group,
    const label size,
    const label start
)
:
    GeoMesh<polyMesh>(mesh.mesh()),
    mesh_(mesh),
    group_(group),
    size_(size),
    start_(start),
    index_(mesh.subMeshIndex())
{}


Foam::LagrangianSubMesh::LagrangianSubMesh
(
    const LagrangianMesh& mesh,
    const labelList& groupOffsets,
    const LagrangianGroup group
)
:
    GeoMesh<polyMesh>(mesh.mesh()),
    mesh_(mesh),
    group_(group),
    size_
    (
        groupOffsets[static_cast<label>(group) + 1]
      - groupOffsets[static_cast<label>(group)]
    ),
    start_(groupOffsets[static_cast<label>(group)]),
    index_(mesh.subMeshIndex())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianSubMesh::~LagrangianSubMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::vectorField> Foam::LagrangianSubMesh::nf
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    tmp<vectorField> tnf(new vectorField(subMesh.size()));
    vectorField& nf = tnf.ref();

    forAll(subMesh, subi)
    {
        nf[subi] =
            tracking::faceNormalAndDisplacement
            (
                subMesh.mesh().mesh(),
                subMesh.mesh().coordinates()[subi + subMesh.start()],
                subMesh.mesh().celli()[subi + subMesh.start()],
                subMesh.mesh().facei()[subi + subMesh.start()],
                subMesh.mesh().faceTrii()[subi + subMesh.start()],
                fraction[subi]
            ).first();
    }

    return tnf;
}


template<>
Foam::tmp<Foam::LagrangianSubVectorField> Foam::LagrangianSubMesh::nf
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    return
        LagrangianSubVectorField::New
        (
            "nf:" + Foam::name(subMesh.group()),
            subMesh,
            dimless,
            nf<vectorField>(fraction)
        );
}


template<>
Foam::tmp<Foam::vectorField> Foam::LagrangianSubMesh::nf
(
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    tmp<vectorField> tnf(new vectorField(size()));
    vectorField& nf = tnf.ref();

    forAll(*this, subi)
    {
        nf[subi] =
            tracking::faceNormalAndDisplacement
            (
                mesh().mesh(),
                mesh().coordinates()[subi + start()],
                mesh().celli()[subi + start()],
                mesh().facei()[subi + start()],
                mesh().faceTrii()[subi + start()],
                fraction[subi + start()]
            ).first();
    }

    return tnf;
}


template<>
Foam::tmp<Foam::LagrangianSubVectorField> Foam::LagrangianSubMesh::nf
(
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    return
        LagrangianSubVectorField::New
        (
            "nf:" + Foam::name(group()),
            *this,
            dimless,
            nf<vectorField>(fraction)
        );
}


template<>
Foam::tmp<Foam::vectorField> Foam::LagrangianSubMesh::Uf
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    tmp<vectorField> tUf(new vectorField(subMesh.size()));
    vectorField& Uf = tUf.ref();

    forAll(subMesh, subi)
    {
        Uf[subi] =
            tracking::faceNormalAndDisplacement
            (
                subMesh.mesh().mesh(),
                subMesh.mesh().coordinates()[subi + subMesh.start()],
                subMesh.mesh().celli()[subi + subMesh.start()],
                subMesh.mesh().facei()[subi + subMesh.start()],
                subMesh.mesh().faceTrii()[subi + subMesh.start()],
                fraction[subi]
            ).second()/subMesh.mesh().time().deltaTValue();
    }

    return tUf;
}


template<>
Foam::tmp<Foam::LagrangianSubVectorField> Foam::LagrangianSubMesh::Uf
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    return
        LagrangianSubVectorField::New
        (
            "Uf:" + Foam::name(subMesh.group()),
            subMesh,
            dimVelocity,
            Uf<vectorField>(fraction)
        );
}


template<>
Foam::tmp<Foam::vectorField> Foam::LagrangianSubMesh::Uf
(
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    tmp<vectorField> tUf(new vectorField(size()));
    vectorField& Uf = tUf.ref();

    forAll(*this, subi)
    {
        Uf[subi] =
            tracking::faceNormalAndDisplacement
            (
                mesh().mesh(),
                mesh().coordinates()[subi + start()],
                mesh().celli()[subi + start()],
                mesh().facei()[subi + start()],
                mesh().faceTrii()[subi + start()],
                fraction[subi + start()]
            ).second()/mesh().time().deltaTValue();
    }

    return tUf;
}


template<>
Foam::tmp<Foam::LagrangianSubVectorField> Foam::LagrangianSubMesh::Uf
(
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    return
        LagrangianSubVectorField::New
        (
            "Uf:" + Foam::name(group()),
            *this,
            dimVelocity,
            Uf<vectorField>(fraction)
        );
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::LagrangianSubMesh::operator+=(const LagrangianSubMesh& subMesh)
{
    if (&mesh_ != &subMesh.mesh_)
    {
        FatalErrorInFunction
            << "Cannot combine sub-meshes which relate to different meshes"
            << exit(FatalError);
    }

    if (group_ != subMesh.group_)
    {
        FatalErrorInFunction
            << "Cannot combine sub-meshes with different groups "
            << exit(FatalError);
    }

    if (size_ + start_ != subMesh.start_)
    {
        FatalErrorInFunction
            << "Cannot combine sub-meshes that are not contiguous"
            << exit(FatalError);
    }

    size_ += subMesh.size_;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "sampledSurface.H"
#include "polyMesh.H"
#include "demandDrivenData.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurface, 0);
    defineRunTimeSelectionTable(sampledSurface, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSurface::makeSf() const
{
    // It is an error to recalculate if the pointer is already set
    if (SfPtr_)
    {
        FatalErrorInFunction
            << "face area vectors already exist"
            << abort(FatalError);
    }

    const faceList& theFaces = faces();
    SfPtr_ = new vectorField(theFaces.size());

    vectorField& values = *SfPtr_;
    forAll(theFaces, facei)
    {
        values[facei] = theFaces[facei].area(points());
    }
}


void Foam::sampledSurface::makeMagSf() const
{
    // It is an error to recalculate if the pointer is already set
    if (magSfPtr_)
    {
        FatalErrorInFunction
            << "mag face areas already exist"
            << abort(FatalError);
    }

    const faceList& theFaces = faces();
    magSfPtr_ = new scalarField(theFaces.size());

    scalarField& values = *magSfPtr_;
    forAll(theFaces, facei)
    {
        values[facei] = theFaces[facei].mag(points());
    }
}


void Foam::sampledSurface::makeCf() const
{
    // It is an error to recalculate if the pointer is already set
    if (CfPtr_)
    {
        FatalErrorInFunction
            << "face centres already exist"
            << abort(FatalError);
    }

    const faceList& theFaces = faces();
    CfPtr_ = new vectorField(theFaces.size());

    vectorField& values = *CfPtr_;
    forAll(theFaces, facei)
    {
        values[facei] = theFaces[facei].centre(points());
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sampledSurface::clearGeom() const
{
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CfPtr_);
    area_ = -1;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sampledSurface> Foam::sampledSurface::New
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const word sampleType(dict.lookup("type"));

    if (debug)
    {
        Info<< "Selecting sampledType " << sampleType << endl;
    }

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(sampleType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown sample type "
            << sampleType << nl << nl
            << "Valid sample types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<sampledSurface>(cstrIter()(name, mesh, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurface::sampledSurface
(
    const word& name,
    const polyMesh& mesh,
    const bool interpolate
)
:
    name_(name),
    mesh_(mesh),
    interpolate_(interpolate),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CfPtr_(nullptr),
    area_(-1)
{}


Foam::sampledSurface::sampledSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    interpolate_(dict.lookupOrDefault("interpolate", false)),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CfPtr_(nullptr),
    area_(-1)
{
    dict.readIfPresent("name", name_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurface::~sampledSurface()
{
    clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::sampledSurface::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}


const Foam::scalarField& Foam::sampledSurface::magSf() const
{
    if (!magSfPtr_)
    {
        makeMagSf();
    }

    return *magSfPtr_;
}


const Foam::vectorField& Foam::sampledSurface::Cf() const
{
    if (!CfPtr_)
    {
        makeCf();
    }

    return *CfPtr_;
}


Foam::scalar Foam::sampledSurface::area() const
{
    if (area_ < 0)
    {
        area_ = sum(magSf());
        reduce(area_, sumOp<scalar>());
    }

    return area_;
}


Foam::tmp<Foam::scalarField> Foam::sampledSurface::sample
(
    const surfaceScalarField& sField
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::vectorField> Foam::sampledSurface::sample
(
    const surfaceVectorField& sField
) const
{
    NotImplemented;
    return tmp<vectorField>(nullptr);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurface::sample
(
    const surfaceSphericalTensorField& sField
) const
{
    NotImplemented;
    return tmp<sphericalTensorField>(nullptr);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurface::sample
(
    const surfaceSymmTensorField& sField
) const
{
    NotImplemented;
    return tmp<symmTensorField>(nullptr);
}


Foam::tmp<Foam::tensorField> Foam::sampledSurface::sample
(
    const surfaceTensorField& sField
) const
{
    NotImplemented;
    return tmp<tensorField>(nullptr);
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sampledSurface::project(const Field<scalar>& field) const
{
    tmp<Field<scalar>> tRes(new Field<scalar>(faces().size()));
    Field<scalar>& res = tRes.ref();

    forAll(faces(), facei)
    {
        res[facei] = field[facei];
    }

    return tRes;
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sampledSurface::project(const Field<vector>& field) const
{
    tmp<Field<scalar>> tRes(new Field<scalar>(faces().size()));
    project(tRes.ref(), field);
    return tRes;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledSurface::project(const Field<sphericalTensor>& field) const
{
    tmp<Field<vector>> tRes(new Field<vector>(faces().size()));
    project(tRes.ref(), field);
    return tRes;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledSurface::project(const Field<symmTensor>& field) const
{
    tmp<Field<vector>> tRes(new Field<vector>(faces().size()));
    project(tRes.ref(), field);
    return tRes;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledSurface::project(const Field<tensor>& field) const
{
    tmp<Field<vector>> tRes(new Field<vector>(faces().size()));
    project(tRes.ref(), field);
    return tRes;
}


void Foam::sampledSurface::print(Ostream& os) const
{
    os << type();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream &os, const sampledSurface& s)
{
    s.print(os);
    os.check("Ostream& operator<<(Ostream&, const sampledSurface&");
    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "cylindrical.H"
#include "axesRotation.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "tensorIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylindrical, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        cylindrical,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        cylindrical,
        objectRegistry
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tensor Foam::cylindrical::R(const vector& dir) const
{
    const vector e3 = e3_/mag(e3_);
    const vector r = dir - (dir & e3)*e3;

    if (mag(r) < small)
    {
        // If the cell centre is on the axis choose any radial direction
        return axesRotation(e3, perpendicular(e3)).R();
    }
    else
    {
        return axesRotation(e3, dir).R();
    }
}


void Foam::cylindrical::init
(
    const objectRegistry& obr,
    const List<label>& cells
)
{
    const polyMesh& mesh = refCast<const polyMesh>(obr);

    Rptr_.reset(new tensorField(cells.size()));

    updateCells(mesh, cells);
}


void Foam::cylindrical::init(const objectRegistry& obr)
{
    const polyMesh& mesh = refCast<const polyMesh>(obr);

    Rptr_.reset(new tensorField(mesh.nCells()));

    const vectorField& cc = mesh.cellCentres();

    tensorField& R = Rptr_();
    forAll(cc, celli)
    {
        vector dir = cc[celli] - origin_;
        dir /= mag(dir) + vSmall;

        R[celli] = this->R(dir);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindrical::cylindrical
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    Rptr_(),
    origin_(point::zero),
    e3_(Zero)
{
    // If origin is specified in the coordinateSystem
    if (dict.parent().found("origin"))
    {
        dict.parent().lookup("origin") >> origin_;
    }

    // rotation axis
    dict.lookup("e3") >> e3_;

    init(obr);
}


Foam::cylindrical::cylindrical
(
    const objectRegistry& obr,
    const vector& axis,
    const point& origin
)
:
    Rptr_(),
    origin_(origin),
    e3_(axis)
{
    init(obr);
}


Foam::cylindrical::cylindrical
(
    const objectRegistry& obr,
    const vector& axis,
    const point& origin,
    const List<label>& cells
)
:
    Rptr_(),
    origin_(origin),
    e3_(axis)
{
    init(obr, cells);
}


Foam::cylindrical::cylindrical(const dictionary& dict)
:
    Rptr_(),
    origin_(),
    e3_()
{
    FatalErrorInFunction
        << " cylindrical can not be constructed from dictionary "
        << " use the constructor : "
           "("
           "    const dictionary&, const objectRegistry&"
           ")"
        << exit(FatalIOError);
}


Foam::cylindrical::cylindrical(const tensorField& R)
:
    Rptr_(),
    origin_(Zero),
    e3_(Zero)
{
    Rptr_() = R;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cylindrical::clear()
{
    if (!Rptr_.empty())
    {
        Rptr_.clear();
    }
}


void Foam::cylindrical::updateCells
(
    const polyMesh& mesh,
    const labelList& cells
)
{
    const vectorField& cc = mesh.cellCentres();
    tensorField& R = Rptr_();

    forAll(cells, i)
    {
        label celli = cells[i];
        vector dir = cc[celli] - origin_;
        dir /= mag(dir) + vSmall;

        R[i] = this->R(dir);
    }
}


Foam::tmp<Foam::vectorField> Foam::cylindrical::transform
(
    const vectorField& vf
) const
{
    if (Rptr_->size() != vf.size())
    {
        FatalErrorInFunction
            << "vectorField st has different size to tensorField "
            << abort(FatalError);
    }

    return (Rptr_() & vf);
}


Foam::vector Foam::cylindrical::transform(const vector& v) const
{
    NotImplemented;
    return Zero;
}


Foam::vector Foam::cylindrical::transform
(
    const vector& v,
    const label cmptI
) const
{
    return (Rptr_()[cmptI] & v);
}


Foam::tmp<Foam::vectorField> Foam::cylindrical::invTransform
(
    const vectorField& vf
) const
{
    return (Rptr_().T() & vf);
}


Foam::vector Foam::cylindrical::invTransform(const vector& v) const
{
    NotImplemented;
    return Zero;
}


Foam::vector Foam::cylindrical::invTransform
(
    const vector& v,
    const label cmptI
) const
{
    return (Rptr_()[cmptI].T() & v);
}


Foam::tmp<Foam::tensorField> Foam::cylindrical::transformTensor
(
    const tensorField& tf
) const
{
    if (Rptr_->size() != tf.size())
    {
        FatalErrorInFunction
            << "tensorField st has different size to tensorField Tr"
            << abort(FatalError);
    }
    return (Rptr_() & tf & Rptr_().T());
}


Foam::tensor Foam::cylindrical::transformTensor
(
    const tensor& t
) const
{
    NotImplemented;

    return Zero;
}


Foam::tmp<Foam::tensorField> Foam::cylindrical::transformTensor
(
    const tensorField& tf,
    const labelList& cellMap
) const
{
    if (cellMap.size() != tf.size())
    {
        FatalErrorInFunction
            << "tensorField tf has different size to tensorField Tr"
            << abort(FatalError);
    }

    const tensorField& R = Rptr_();
    const tensorField Rtr(R.T());
    tmp<tensorField> tt(new tensorField(cellMap.size()));
    tensorField& t = tt.ref();
    forAll(cellMap, i)
    {
        t[i] = R[i] & tf[i] & Rtr[i];
    }

    return tt;
}


Foam::tmp<Foam::symmTensorField> Foam::cylindrical::transformVector
(
    const vectorField& vf
) const
{
    if (Rptr_->size() != vf.size())
    {
        FatalErrorInFunction
            << "tensorField vf has different size to tensorField Tr"
            << abort(FatalError);
    }

    tmp<symmTensorField> tfld(new symmTensorField(Rptr_->size()));
    symmTensorField& fld = tfld.ref();

    const tensorField& R = Rptr_();
    forAll(fld, i)
    {
        fld[i] = transformPrincipal(R[i], vf[i]);
    }
    return tfld;
}


Foam::symmTensor Foam::cylindrical::transformVector
(
    const vector& v
) const
{
    NotImplemented;
    return Zero;
}


void Foam::cylindrical::write(Ostream& os) const
{
     writeEntry(os, "e3", e3());
}


// ************************************************************************* //

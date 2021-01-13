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

#include "cylindrical.H"
#include "axesRotation.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylindrical, 0);
    addToRunTimeSelectionTable(coordinateRotation, cylindrical, dictionary);
    addToRunTimeSelectionTable(coordinateRotation, cylindrical, points);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tensor Foam::cylindrical::R(const vector& p) const
{
    vector dir = p - origin_;
    dir /= mag(dir) + vSmall;

    const vector axis = axis_/mag(axis_);
    const vector r = dir - (dir & axis)*axis;

    if (mag(r) < small)
    {
        // If the point is on the axis choose any radial direction
        return axesRotation(axis, perpendicular(axis)).R();
    }
    else
    {
        return axesRotation(axis, dir).R();
    }

    return tensor(r, axis^r, axis);
}


void Foam::cylindrical::init(const UList<vector>& points)
{
    Rptr_.reset(new tensorField(points.size()));
    tensorField& R = Rptr_();

    forAll(points, i)
    {
        R[i] = this->R(points[i]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindrical::cylindrical
(
    const vector& axis,
    const point& origin,
    const UList<vector>& points
)
:
    Rptr_(),
    origin_(origin),
    axis_(axis)
{
    init(points);
}


Foam::cylindrical::cylindrical(const dictionary& dict)
:
    Rptr_(),
    origin_
    (
        dict.parent().found("origin")
      ? dict.parent().lookup("origin")
      : dict.lookup("origin")
    ),
    axis_(dict.lookupBackwardsCompatible<vector>({"axis", "e3"}))
{}


Foam::cylindrical::cylindrical
(
    const dictionary& dict,
    const UList<vector>& points
)
:
    cylindrical(dict)
{
    init(points);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cylindrical::updatePoints(const UList<vector>& points)
{
    if (Rptr_.valid())
    {
        Rptr_().setSize(points.size());
    }
    else
    {
        Rptr_.reset(new tensorField(points.size()));
    }

    tensorField& R = Rptr_();

    forAll(points, i)
    {
        R[i] = this->R(points[i]);
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
    if (Rptr_->size() != vf.size())
    {
        FatalErrorInFunction
            << "vectorField st has different size to tensorField "
            << abort(FatalError);
    }

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
     writeEntry(os, "axis", axis());
}


// ************************************************************************* //

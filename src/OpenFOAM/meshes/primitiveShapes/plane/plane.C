/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "plane.H"
#include "tensor.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plane::calcPntAndVec
(
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar d
)
{
    normal_ = vector(a, b, c);

    const scalar magNormal = mag(normal_);

    // Normalise the normal if possible. Set to invalid if not.
    if (magNormal > 0)
    {
        normal_ /= magNormal;
    }
    else
    {
        normal_ = vector::zero;
    }

    // Construct the point if possible. Set to far away if not.
    if (magNormal > mag(d)*vSmall)
    {
        point_ = - d/magNormal*normal_;
    }
    else
    {
        point_ = point::max;
    }
}


void Foam::plane::calcPntAndVec
(
    const point& point1,
    const point& point2,
    const point& point3
)
{
    normal_ = normalised((point1 - point2) ^ (point2 - point3));

    point_ = (point1 + point2 + point3)/3;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plane::plane(const vector& normalVector)
:
    normal_(normalised(normalVector)),
    point_(Zero)
{}


Foam::plane::plane(const point& basePoint, const vector& normalVector)
:
    normal_(normalised(normalVector)),
    point_(basePoint)
{}


Foam::plane::plane
(
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar d
)
{
    calcPntAndVec(a, b, c, d);
}


Foam::plane::plane
(
    const point& point1,
    const point& point2,
    const point& point3
)
{
    calcPntAndVec(point1, point2, point3);
}


Foam::plane::plane(const dictionary& dict)
:
    normal_(Zero),
    point_(Zero)
{
    const word planeType(dict.lookup("planeType"));

    const dictionary& subDict = dict.optionalSubDict(planeType + "Dict");

    if (planeType == "planeEquation")
    {
        calcPntAndVec
        (
            subDict.lookup<scalar>("a"),
            subDict.lookup<scalar>("b"),
            subDict.lookup<scalar>("c"),
            subDict.lookup<scalar>("d")
        );
    }
    else if (planeType == "embeddedPoints")
    {
        calcPntAndVec
        (
            subDict.lookup<point>("point1"),
            subDict.lookup<point>("point2"),
            subDict.lookup<point>("point3")
        );
    }
    else if (planeType == "pointAndNormal")
    {
        point_ =
            subDict.lookupBackwardsCompatible<point>
            (
                {"point", "basePoint"}
            );

        normal_ =
            normalised
            (
                subDict.lookupBackwardsCompatible<vector>
                (
                    {"normal", "normalVector"}
                )
            );
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Invalid plane type: " << planeType << nl
            << "Valid options include: "
            << "planeEquation, embeddedPoints and pointAndNormal"
            << abort(FatalIOError);
    }

    if (normal_ == vector::zero)
    {
        FatalIOErrorInFunction(subDict)
            << "Plane normal has zero length"
            << exit(FatalIOError);
    }

    if (point_ == point::max)
    {
        FatalIOErrorInFunction(subDict)
            << "Plane is too far from the origin"
            << exit(FatalIOError);
    }
}


Foam::plane::plane(Istream& is)
:
    normal_(normalised(vector(is))),
    point_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::scalar, 4> Foam::plane::planeCoeffs() const
{
    FixedList<scalar, 4> C(4);

    scalar magX = mag(normal_.x());
    scalar magY = mag(normal_.y());
    scalar magZ = mag(normal_.z());

    if (magX > magY)
    {
        if (magX > magZ)
        {
            C[0] = 1;
            C[1] = normal_.y()/normal_.x();
            C[2] = normal_.z()/normal_.x();
        }
        else
        {
            C[0] = normal_.x()/normal_.z();
            C[1] = normal_.y()/normal_.z();
            C[2] = 1;
        }
    }
    else
    {
        if (magY > magZ)
        {
            C[0] = normal_.x()/normal_.y();
            C[1] = 1;
            C[2] = normal_.z()/normal_.y();
        }
        else
        {
            C[0] = normal_.x()/normal_.z();
            C[1] = normal_.y()/normal_.z();
            C[2] = 1;
        }
    }

    C[3] = - C[0] * point_.x()
           - C[1] * point_.y()
           - C[2] * point_.z();

    return C;
}


Foam::point Foam::plane::aPoint() const
{
    // Perturb base point
    const point& refPt = refPoint();

    // ax + by + cz + d = 0
    FixedList<scalar, 4> planeCoeffs = this->planeCoeffs();

    const scalar perturbX = refPt.x() + 1e-3;
    const scalar perturbY = refPt.y() + 1e-3;
    const scalar perturbZ = refPt.z() + 1e-3;

    if (mag(planeCoeffs[2]) < small)
    {
        if (mag(planeCoeffs[1]) < small)
        {
            const scalar x =
                -1.0
                *(
                     planeCoeffs[3]
                   + planeCoeffs[1]*perturbY
                   + planeCoeffs[2]*perturbZ
                 )/planeCoeffs[0];

            return point
            (
                x,
                perturbY,
                perturbZ
            );
        }

        const scalar y =
            -1.0
            *(
                 planeCoeffs[3]
               + planeCoeffs[0]*perturbX
               + planeCoeffs[2]*perturbZ
             )/planeCoeffs[1];

        return point
        (
            perturbX,
            y,
            perturbZ
        );
    }
    else
    {
        const scalar z =
            -1.0
            *(
                 planeCoeffs[3]
               + planeCoeffs[0]*perturbX
               + planeCoeffs[1]*perturbY
             )/planeCoeffs[2];

        return point
        (
            perturbX,
            perturbY,
            z
        );
    }
}


Foam::point Foam::plane::nearestPoint(const point& p) const
{
    return p - normal_*signedDistance(p);
}


Foam::scalar Foam::plane::distance(const point& p) const
{
    return mag(signedDistance(p));
}


Foam::scalar Foam::plane::signedDistance(const point& p) const
{
    return (p - point_) & normal_;
}


Foam::scalar Foam::plane::normalIntersect
(
    const point& pnt0,
    const vector& dir
) const
{
    const scalar num = (point_ - pnt0) & normal_;
    const scalar den = dir & normal_;

    return mag(den) > mag(num)*vSmall ? num/den : vGreat;
}


Foam::plane::ray Foam::plane::planeIntersect(const plane& plane2) const
{
    // Mathworld plane-plane intersection. Assume there is a point on the
    // intersection line with z=0 and solve the two plane equations
    // for that (now 2x2 equation in x and y)
    // Better: use either z=0 or x=0 or y=0.

    const vector& n1 = normal();
    const vector& n2 = plane2.normal();

    const point& p1 = refPoint();
    const point& p2 = plane2.refPoint();

    scalar n1p1 = n1&p1;
    scalar n2p2 = n2&p2;

    vector dir = n1 ^ n2;

    // Determine zeroed out direction (can be x,y or z) by looking at which
    // has the largest component in dir.
    scalar magX = mag(dir.x());
    scalar magY = mag(dir.y());
    scalar magZ = mag(dir.z());

    direction iZero, i1, i2;

    if (magX > magY)
    {
        if (magX > magZ)
        {
            iZero = 0;
            i1 = 1;
            i2 = 2;
        }
        else
        {
            iZero = 2;
            i1 = 0;
            i2 = 1;
        }
    }
    else
    {
        if (magY > magZ)
        {
            iZero = 1;
            i1 = 2;
            i2 = 0;
        }
        else
        {
            iZero = 2;
            i1 = 0;
            i2 = 1;
        }
    }

    vector pt;

    pt[iZero] = 0;
    pt[i1] = (n2[i2]*n1p1 - n1[i2]*n2p2) / (n1[i1]*n2[i2] - n2[i1]*n1[i2]);
    pt[i2] = (n2[i1]*n1p1 - n1[i1]*n2p2) / (n1[i2]*n2[i1] - n1[i1]*n2[i2]);

    return ray(pt, dir);
}


Foam::point Foam::plane::planePlaneIntersect
(
    const plane& plane2,
    const plane& plane3
) const
{
    FixedList<scalar, 4> coeffs1(planeCoeffs());
    FixedList<scalar, 4> coeffs2(plane2.planeCoeffs());
    FixedList<scalar, 4> coeffs3(plane3.planeCoeffs());

    tensor a
    (
        coeffs1[0],coeffs1[1],coeffs1[2],
        coeffs2[0],coeffs2[1],coeffs2[2],
        coeffs3[0],coeffs3[1],coeffs3[2]
    );

    vector b(coeffs1[3],coeffs2[3],coeffs3[3]);

    return (inv(a) & (-b));
}


Foam::plane::side Foam::plane::sideOfPlane(const point& p) const
{
    const scalar angle((p - point_) & normal_);

    return (angle < 0 ? FLIP : NORMAL);
}


Foam::point Foam::plane::mirror(const point& p) const
{
    const vector mirroredPtDir = p - nearestPoint(p);

    if ((normal() & mirroredPtDir) > 0)
    {
        return p - 2.0*distance(p)*normal();
    }
    else
    {
        return p + 2.0*distance(p)*normal();
    }
}


void Foam::plane::writeDict(Ostream& os) const
{
    writeEntry(os, "planeType", "pointAndNormal");
    os  << indent << "pointAndNormalDict" << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;
    writeEntry(os, "point", point_);
    writeEntry(os, "normal", normal_);
    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const plane& a, const plane& b)
{
    if (a.point_ == b.point_ && a.normal_ == b.normal_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::operator!=(const plane& a, const plane& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const plane& a)
{
    os  << a.normal_ << token::SPACE << a.point_;

    return os;
}


// ************************************************************************* //

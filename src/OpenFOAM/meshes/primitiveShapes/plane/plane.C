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

#include "plane.H"
#include "tensor.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plane::calcPntAndVec(const scalarList& C)
{
    normal_ = vector(C[0], C[1], C[2]);

    const scalar magNormal = mag(normal_);

    if (magNormal == 0)
    {
        FatalErrorInFunction
            << "Plane normal has zero length"
            << abort(FatalError);
    }

    normal_ /= magNormal;

    if (magNormal < mag(C[3])*vSmall)
    {
        FatalErrorInFunction
            << "Plane is too far from the origin"
            << abort(FatalError);
    }

    point_ = - C[3]/magNormal*normal_;
}


void Foam::plane::calcPntAndVec
(
    const point& point1,
    const point& point2,
    const point& point3
)
{
    point_ = (point1 + point2 + point3)/3;
    vector line12 = point1 - point2;
    vector line23 = point2 - point3;

    if
    (
        mag(line12) < vSmall
     || mag(line23) < vSmall
     || mag(point3-point1) < vSmall
    )
    {
        FatalErrorInFunction
            << "Bad points:" << point1 << ' ' << point2 << ' ' << point3
            << abort(FatalError);
    }

    normal_ = line12 ^ line23;
    scalar magUnitVector(mag(normal_));

    if (magUnitVector < vSmall)
    {
        FatalErrorInFunction
            << "Plane normal defined with zero length" << nl
            << "Bad points:" << point1 << ' ' << point2 << ' ' << point3
            << abort(FatalError);
    }

    normal_ /= magUnitVector;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plane::plane(const vector& normalVector)
:
    normal_(normalVector),
    point_(Zero)
{
    scalar magUnitVector(mag(normal_));

    if (magUnitVector > vSmall)
    {
        normal_ /= magUnitVector;
    }
    else
    {
        FatalErrorInFunction
            << "plane normal has zero length. basePoint:" << point_
            << abort(FatalError);
    }
}


Foam::plane::plane(const point& basePoint, const vector& normalVector)
:
    normal_(normalVector),
    point_(basePoint)
{
    scalar magUnitVector(mag(normal_));

    if (magUnitVector > vSmall)
    {
        normal_ /= magUnitVector;
    }
    else
    {
        FatalErrorInFunction
            << "plane normal has zero length. basePoint:" << point_
            << abort(FatalError);
    }
}


Foam::plane::plane(const scalarList& C)
{
    calcPntAndVec(C);
}


Foam::plane::plane
(
    const point& a,
    const point& b,
    const point& c
)
{
    calcPntAndVec(a, b, c);
}


Foam::plane::plane(const dictionary& dict)
:
    normal_(Zero),
    point_(Zero)
{
    const word planeType(dict.lookup("planeType"));

    if (planeType == "planeEquation")
    {
        const dictionary& subDict = dict.subDict("planeEquationDict");
        scalarList C(4);

        C[0] = subDict.lookup<scalar>("a");
        C[1] = subDict.lookup<scalar>("b");
        C[2] = subDict.lookup<scalar>("c");
        C[3] = subDict.lookup<scalar>("d");

        calcPntAndVec(C);

    }
    else if (planeType == "embeddedPoints")
    {
        const dictionary& subDict = dict.subDict("embeddedPointsDict");

        point point1(subDict.lookup("point1"));
        point point2(subDict.lookup("point2"));
        point point3(subDict.lookup("point3"));

        calcPntAndVec(point1, point2, point3);
    }
    else if (planeType == "pointAndNormal")
    {
        const dictionary& subDict = dict.subDict("pointAndNormalDict");

        point_ =
            subDict.found("basePoint")
          ? subDict.lookup("basePoint")
          : subDict.lookup("point");

        normal_ =
            subDict.found("normalVector")
          ? subDict.lookup("normalVector")
          : subDict.lookup("normal");

        normal_ /= mag(normal_);
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Invalid plane type: " << planeType << nl
            << "Valid options include: planeEquation, embeddedPoints and "
            << "pointAndNormal"
            << abort(FatalIOError);
    }
}


Foam::plane::plane(Istream& is)
:
    normal_(is),
    point_(is)
{
    scalar magUnitVector(mag(normal_));

    if (magUnitVector > vSmall)
    {
        normal_ /= magUnitVector;
    }
    else
    {
        FatalErrorInFunction
            << "plane normal has zero length. basePoint:" << point_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vector& Foam::plane::normal() const
{
    return normal_;
}


const Foam::point& Foam::plane::refPoint() const
{
    return point_;
}


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
    return p - normal_*((p - point_) & normal_);
}


Foam::scalar Foam::plane::distance(const point& p) const
{
    return mag((p - point_) & normal_);
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

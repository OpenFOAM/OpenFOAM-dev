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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::quaternion::quaternion()
{}


inline Foam::quaternion::quaternion(const scalar w, const vector& v)
:
    w_(w),
    v_(v)
{}


inline Foam::quaternion::quaternion(const vector& d, const scalar theta)
:
    w_(cos(0.5*theta)),
    v_((sin(0.5*theta)/mag(d))*d)
{}


inline Foam::quaternion::quaternion
(
    const vector& d,
    const scalar cosTheta,
    const bool normalized
)
{
    scalar cosHalfTheta2 = 0.5*(cosTheta + 1);
    w_ = sqrt(cosHalfTheta2);

    if (normalized)
    {
        v_ = sqrt(1 - cosHalfTheta2)*d;
    }
    else
    {
        v_ = (sqrt(1 - cosHalfTheta2)/mag(d))*d;
    }
}


inline Foam::quaternion::quaternion(const scalar w)
:
    w_(w),
    v_(Zero)
{}


inline Foam::quaternion::quaternion(const vector& v)
:
    w_(0),
    v_(v)
{}


inline Foam::quaternion Foam::quaternion::unit(const vector& v)
{
    return quaternion(sqrt(1 - magSqr(v)), v);
}


inline Foam::quaternion::quaternion
(
    const rotationSequence rs,
    const vector& angles
)
{
    switch(rs)
    {
        case ZYX:
            operator=(quaternion(vector(0, 0, 1), angles.x()));
            operator*=(quaternion(vector(0, 1, 0), angles.y()));
            operator*=(quaternion(vector(1, 0, 0), angles.z()));
            break;

        case ZYZ:
            operator=(quaternion(vector(0, 0, 1), angles.x()));
            operator*=(quaternion(vector(0, 1, 0), angles.y()));
            operator*=(quaternion(vector(0, 0, 1), angles.z()));
            break;

        case ZXY:
            operator=(quaternion(vector(0, 0, 1), angles.x()));
            operator*=(quaternion(vector(1, 0, 0), angles.y()));
            operator*=(quaternion(vector(0, 1, 0), angles.z()));
            break;

        case ZXZ:
            operator=(quaternion(vector(0, 0, 1), angles.x()));
            operator*=(quaternion(vector(1, 0, 0), angles.y()));
            operator*=(quaternion(vector(0, 0, 1), angles.z()));
            break;

        case YXZ:
            operator=(quaternion(vector(0, 1, 0), angles.x()));
            operator*=(quaternion(vector(1, 0, 0), angles.y()));
            operator*=(quaternion(vector(0, 0, 1), angles.z()));
            break;

        case YXY:
            operator=(quaternion(vector(0, 1, 0), angles.x()));
            operator*=(quaternion(vector(1, 0, 0), angles.y()));
            operator*=(quaternion(vector(0, 1, 0), angles.z()));
            break;

        case YZX:
            operator=(quaternion(vector(0, 1, 0), angles.x()));
            operator*=(quaternion(vector(0, 0, 1), angles.y()));
            operator*=(quaternion(vector(1, 0, 0), angles.z()));
            break;

        case YZY:
            operator=(quaternion(vector(0, 1, 0), angles.x()));
            operator*=(quaternion(vector(0, 0, 1), angles.y()));
            operator*=(quaternion(vector(0, 1, 0), angles.z()));
            break;

        case XYZ:
            operator=(quaternion(vector(1, 0, 0), angles.x()));
            operator*=(quaternion(vector(0, 1, 0), angles.y()));
            operator*=(quaternion(vector(0, 0, 1), angles.z()));
            break;

        case XYX:
            operator=(quaternion(vector(1, 0, 0), angles.x()));
            operator*=(quaternion(vector(0, 1, 0), angles.y()));
            operator*=(quaternion(vector(1, 0, 0), angles.z()));
            break;

        case XZY:
            operator=(quaternion(vector(1, 0, 0), angles.x()));
            operator*=(quaternion(vector(0, 0, 1), angles.y()));
            operator*=(quaternion(vector(0, 1, 0), angles.z()));
            break;

        case XZX:
            operator=(quaternion(vector(1, 0, 0), angles.x()));
            operator*=(quaternion(vector(0, 0, 1), angles.y()));
            operator*=(quaternion(vector(1, 0, 0), angles.z()));
            break;

        default:
            FatalErrorInFunction
                << "Unknown rotation sequence " << rs << abort(FatalError);
            break;
    }
}


inline Foam::quaternion::quaternion
(
    const tensor& rotationTensor
)
{
    scalar trace =
        rotationTensor.xx()
      + rotationTensor.yy()
      + rotationTensor.zz();

    if (trace > 0)
    {
        scalar s = 0.5/Foam::sqrt(trace + 1.0);

        w_ = 0.25/s;
        v_[0] = (rotationTensor.zy() - rotationTensor.yz())*s;
        v_[1] = (rotationTensor.xz() - rotationTensor.zx())*s;
        v_[2] = (rotationTensor.yx() - rotationTensor.xy())*s;
    }
    else
    {
        if
        (
            rotationTensor.xx() > rotationTensor.yy()
         && rotationTensor.xx() > rotationTensor.zz()
        )
        {
            scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.xx()
              - rotationTensor.yy()
              - rotationTensor.zz()
            );

            w_ = (rotationTensor.zy() - rotationTensor.yz())/s;
            v_[0] = 0.25*s;
            v_[1] = (rotationTensor.xy() + rotationTensor.yx())/s;
            v_[2] = (rotationTensor.xz() + rotationTensor.zx())/s;
        }
        else if
        (
            rotationTensor.yy() > rotationTensor.zz()
        )
        {
            scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.yy()
              - rotationTensor.xx()
              - rotationTensor.zz()
            );

            w_ = (rotationTensor.xz() - rotationTensor.xz())/s;
            v_[0] = (rotationTensor.xy() + rotationTensor.yx())/s;
            v_[1] = 0.25*s;
            v_[2] = (rotationTensor.yz() + rotationTensor.zy())/s;
        }
        else
        {
            scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.zz()
              - rotationTensor.xx()
              - rotationTensor.yy()
            );

            w_ = (rotationTensor.yx() - rotationTensor.xy())/s;
            v_[0] = (rotationTensor.xz() + rotationTensor.zx())/s;
            v_[1] = (rotationTensor.yz() + rotationTensor.zy())/s;
            v_[2] = 0.25*s;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::quaternion::w() const
{
    return w_;
}


inline const Foam::vector& Foam::quaternion::v() const
{
    return v_;
}


inline Foam::scalar& Foam::quaternion::w()
{
    return w_;
}


inline Foam::vector& Foam::quaternion::v()
{
    return v_;
}


inline Foam::quaternion Foam::quaternion::normalized() const
{
    return *this/mag(*this);
}


inline void Foam::quaternion::normalize()
{
    operator/=(mag(*this));
}


inline Foam::quaternion Foam::quaternion::mulq0v(const vector& u) const
{
    return quaternion(-(v() & u), w()*u + (v() ^ u));
}


inline Foam::vector Foam::quaternion::transform(const vector& u) const
{
    return (mulq0v(u)*conjugate(*this)).v();
}


inline Foam::vector Foam::quaternion::invTransform(const vector& u) const
{
    return (conjugate(*this).mulq0v(u)*(*this)).v();
}


inline Foam::quaternion Foam::quaternion::transform(const quaternion& q) const
{
    return Foam::normalize((*this)*q);
}


inline Foam::quaternion Foam::quaternion::invTransform
(
    const quaternion& q
) const
{
    return Foam::normalize(conjugate(*this)*q);
}


inline Foam::tensor Foam::quaternion::R() const
{
    const scalar w2 = sqr(w());
    const scalar x2 = sqr(v().x());
    const scalar y2 = sqr(v().y());
    const scalar z2 = sqr(v().z());

    const scalar txy = 2*v().x()*v().y();
    const scalar twz = 2*w()*v().z();
    const scalar txz = 2*v().x()*v().z();
    const scalar twy = 2*w()*v().y();
    const scalar tyz = 2*v().y()*v().z();
    const scalar twx = 2*w()*v().x();

    return tensor
    (
        w2 + x2 - y2 - z2,  txy - twz,          txz + twy,
        txy + twz,          w2 - x2 + y2 - z2,  tyz - twx,
        txz - twy,          tyz + twx,          w2 - x2 - y2 + z2
    );
}


inline Foam::vector Foam::quaternion::twoAxes
(
    const scalar t11,
    const scalar t12,
    const scalar c2,
    const scalar t31,
    const scalar t32
)
{
    return vector(atan2(t11, t12), acos(c2), atan2(t31, t32));
}


inline Foam::vector Foam::quaternion::threeAxes
(
    const scalar t11,
    const scalar t12,
    const scalar s2,
    const scalar t31,
    const scalar t32
)
{
    return vector(atan2(t11, t12), asin(s2), atan2(t31, t32));
}


inline Foam::vector Foam::quaternion::eulerAngles
(
    const rotationSequence rs
) const
{
    const scalar w2 = sqr(w());
    const scalar x2 = sqr(v().x());
    const scalar y2 = sqr(v().y());
    const scalar z2 = sqr(v().z());

    switch(rs)
    {
        case ZYX:
            return threeAxes
            (
                2*(v().x()*v().y() + w()*v().z()),
                w2 + x2 - y2 - z2,
                2*(w()*v().y() - v().x()*v().z()),
                2*(v().y()*v().z() + w()*v().x()),
                w2 - x2 - y2 + z2
            );
            break;

        case ZYZ:
            return twoAxes
            (
                2*(v().y()*v().z() - w()*v().x()),
                2*(v().x()*v().z() + w()*v().y()),
                w2 - x2 - y2 + z2,
                2*(v().y()*v().z() + w()*v().x()),
                2*(w()*v().y() - v().x()*v().z())
            );
            break;

        case ZXY:
            return threeAxes
            (
                2*(w()*v().z() - v().x()*v().y()),
                w2 - x2 + y2 - z2,
                2*(v().y()*v().z() + w()*v().x()),
                2*(w()*v().y() - v().x()*v().z()),
                w2 - x2 - y2 + z2
            );
            break;

        case ZXZ:
            return twoAxes
            (
                2*(v().x()*v().z() + w()*v().y()),
                2*(w()*v().x() - v().y()*v().z()),
                w2 - x2 - y2 + z2,
                2*(v().x()*v().z() - w()*v().y()),
                2*(v().y()*v().z() + w()*v().x())
            );
            break;

        case YXZ:
            return threeAxes
            (
                2*(v().x()*v().z() + w()*v().y()),
                w2 - x2 - y2 + z2,
                2*(w()*v().x() - v().y()*v().z()),
                2*(v().x()*v().y() + w()*v().z()),
                w2 - x2 + y2 - z2
            );
            break;

        case YXY:
            return twoAxes
            (
                2*(v().x()*v().y() - w()*v().z()),
                2*(v().y()*v().z() + w()*v().x()),
                w2 - x2 + y2 - z2,
                2*(v().x()*v().y() + w()*v().z()),
                2*(w()*v().x() - v().y()*v().z())
            );
            break;

        case YZX:
            return threeAxes
            (
                2*(w()*v().y() - v().x()*v().z()),
                w2 + x2 - y2 - z2,
                2*(v().x()*v().y() + w()*v().z()),
                2*(w()*v().x() - v().y()*v().z()),
                w2 - x2 + y2 - z2
            );
            break;

        case YZY:
            return twoAxes
            (
                2*(v().y()*v().z() + w()*v().x()),
                2*(w()*v().z() - v().x()*v().y()),
                w2 - x2 + y2 - z2,
                2*(v().y()*v().z() - w()*v().x()),
                2*(v().x()*v().y() + w()*v().z())
            );
            break;

        case XYZ:
            return threeAxes
            (
                2*(w()*v().x() - v().y()*v().z()),
                w2 - x2 - y2 + z2,
                2*(v().x()*v().z() + w()*v().y()),
                2*(w()*v().z() - v().x()*v().y()),
                w2 + x2 - y2 - z2
            );
            break;

        case XYX:
            return twoAxes
            (
                2*(v().x()*v().y() + w()*v().z()),
                2*(w()*v().y() - v().x()*v().z()),
                w2 + x2 - y2 - z2,
                2*(v().x()*v().y() - w()*v().z()),
                2*(v().x()*v().z() + w()*v().y())
            );
            break;

        case XZY:
            return threeAxes
            (
                2*(v().y()*v().z() + w()*v().x()),
                w2 - x2 + y2 - z2,
                2*(w()*v().z() - v().x()*v().y()),
                2*(v().x()*v().z() + w()*v().y()),
                w2 + x2 - y2 - z2
            );
            break;

        case XZX:
            return twoAxes
            (
                2*(v().x()*v().z() - w()*v().y()),
                2*(v().x()*v().y() + w()*v().z()),
                w2 + x2 - y2 - z2,
                2*(v().x()*v().z() + w()*v().y()),
                2*(w()*v().z() - v().x()*v().y())
            );
            break;
        default:
            FatalErrorInFunction
                << "Unknown rotation sequence " << rs << abort(FatalError);
            return Zero;
            break;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::quaternion::operator+=(const quaternion& q)
{
    w_ += q.w_;
    v_ += q.v_;
}

inline void Foam::quaternion::operator-=(const quaternion& q)
{
    w_ -= q.w_;
    v_ -= q.v_;
}

inline void Foam::quaternion::operator*=(const quaternion& q)
{
    scalar w0 = w();
    w() = w()*q.w() - (v() & q.v());
    v() = w0*q.v() + q.w()*v() + (v() ^ q.v());
}

inline void Foam::quaternion::operator/=(const quaternion& q)
{
    return operator*=(inv(q));
}


inline void Foam::quaternion::operator=(const scalar s)
{
    w_ = s;
}


inline void Foam::quaternion::operator=(const vector& v)
{
    v_ = v;
}


inline void Foam::quaternion::operator*=(const scalar s)
{
    w_ *= s;
    v_ *= s;
}

inline void Foam::quaternion::operator/=(const scalar s)
{
    w_ /= s;
    v_ /= s;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::magSqr(const quaternion& q)
{
    return magSqr(q.w()) + magSqr(q.v());
}


inline Foam::scalar Foam::mag(const quaternion& q)
{
    return sqrt(magSqr(q));
}


inline Foam::quaternion Foam::conjugate(const quaternion& q)
{
    return quaternion(q.w(), -q.v());
}


inline Foam::quaternion Foam::inv(const quaternion& q)
{
    scalar magSqrq = magSqr(q);
    return quaternion(q.w()/magSqrq, -q.v()/magSqrq);
}


inline Foam::quaternion Foam::normalize(const quaternion& q)
{
    return q/mag(q);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const quaternion& q1, const quaternion& q2)
{
    return (equal(q1.w(), q2.w()) && equal(q1.v(), q2.v()));
}


inline bool Foam::operator!=(const quaternion& q1, const quaternion& q2)
{
    return !operator==(q1, q2);
}


inline Foam::quaternion Foam::operator+
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return quaternion(q1.w() + q2.w(), q1.v() + q2.v());
}


inline Foam::quaternion Foam::operator-(const quaternion& q)
{
    return quaternion(-q.w(), -q.v());
}


inline Foam::quaternion Foam::operator-
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return quaternion(q1.w() - q2.w(), q1.v() - q2.v());
}


inline Foam::scalar Foam::operator&(const quaternion& q1, const quaternion& q2)
{
    return q1.w()*q2.w() + (q1.v() & q2.v());
}


inline Foam::quaternion Foam::operator*
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return quaternion
    (
        q1.w()*q2.w() - (q1.v() & q2.v()),
        q1.w()*q2.v() + q2.w()*q1.v() + (q1.v() ^ q2.v())
    );
}


inline Foam::quaternion Foam::operator/
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return q1*inv(q2);
}


inline Foam::quaternion Foam::operator*(const scalar s, const quaternion& q)
{
    return quaternion(s*q.w(), s*q.v());
}


inline Foam::quaternion Foam::operator*(const quaternion& q, const scalar s)
{
    return quaternion(s*q.w(), s*q.v());
}


inline Foam::quaternion Foam::operator/(const quaternion& q, const scalar s)
{
    return quaternion(q.w()/s, q.v()/s);
}


// ************************************************************************* //

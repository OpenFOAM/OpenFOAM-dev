/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "triad.H"
#include "transform.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::triad::vsType::typeName = "triad";

template<>
const char* const Foam::triad::vsType::componentNames[] = {"x", "y", "z"};

template<>
const Foam::Vector<Foam::vector> Foam::triad::vsType::zero
(
    triad::uniform(vector::uniform(0))
);

template<>
const Foam::Vector<Foam::vector> Foam::triad::vsType::one
(
    triad::uniform(vector::uniform(1))
);

template<>
const Foam::Vector<Foam::vector> Foam::triad::vsType::max
(
    triad::uniform(vector::uniform(vGreat))
);

template<>
const Foam::Vector<Foam::vector> Foam::triad::vsType::min
(
    triad::uniform(vector::uniform(-vGreat))
);

template<>
const Foam::Vector<Foam::vector> Foam::triad::vsType::rootMax
(
    triad::uniform(vector::uniform(rootVGreat))
);

template<>
const Foam::Vector<Foam::vector> Foam::triad::vsType::rootMin
(
    triad::uniform(vector::uniform(-rootVGreat))
);

const Foam::triad Foam::triad::I
(
    vector(1, 0, 0),
    vector(0, 1, 0),
    vector(0, 0, 1)
);

const Foam::triad Foam::triad::unset
(
    triad::uniform(vector::uniform(vGreat))
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triad::triad(const quaternion& q)
{
    tensor Rt(q.R().T());
    x() = Rt.x();
    y() = Rt.y();
    z() = Rt.z();
}


Foam::triad::triad(const tensor& t)
{
    x() = t.x();
    y() = t.y();
    z() = t.z();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triad::orthogonalise()
{
    // Hack for 2D z-slab cases
    // if (!set(2))
    // {
    //     operator[](2) = vector(0, 0, 1);
    // }

    // If only two of the axes are set, set the third
    if (set(0) && set(1) && !set(2))
    {
        operator[](2) = orthogonal(operator[](0), operator[](1));
    }
    else if (set(0) && set(2) && !set(1))
    {
        operator[](1) = orthogonal(operator[](0), operator[](2));
    }
    else if (set(1) && set(2) && !set(0))
    {
        operator[](0) = orthogonal(operator[](1), operator[](2));
    }

    // If all the axes are set
    if (set())
    {
        for (int i=0; i<2; i++)
        {
            const scalar o01
            (
                (set(0) && set(1)) ? mag(operator[](0) & operator[](1)) : vGreat
            );
            const scalar o02
            (
                (set(0) && set(2)) ? mag(operator[](0) & operator[](2)) : vGreat
            );
            const scalar o12
            (
                (set(1) && set(2)) ? mag(operator[](1) & operator[](2)) : vGreat
            );

            if (o01 < o02 && o01 < o12)
            {
                operator[](2) = orthogonal(operator[](0), operator[](1));

                // if (o02 < o12)
                // {
                //     operator[](1) = orthogonal(operator[](0), operator[](2));
                // }
                // else
                // {
                //     operator[](0) = orthogonal(operator[](1), operator[](2));
                // }
            }
            else if (o02 < o12)
            {
                operator[](1) = orthogonal(operator[](0), operator[](2));

                // if (o01 < o12)
                // {
                //     operator[](2) = orthogonal(operator[](0), operator[](1));
                // }
                // else
                // {
                //     operator[](0) = orthogonal(operator[](1), operator[](2));
                // }
            }
            else
            {
                operator[](0) = orthogonal(operator[](1), operator[](2));

                // if (o02 < o01)
                // {
                //     operator[](1) = orthogonal(operator[](0), operator[](2));
                // }
                // else
                // {
                //     operator[](2) = orthogonal(operator[](0), operator[](1));
                // }
            }
        }
    }
}


void Foam::triad::operator+=(const triad& t2)
{
    bool preset[3];

    for (direction i=0; i<3; i++)
    {
        if (t2.set(i) && !set(i))
        {
            operator[](i) = t2.operator[](i);
            preset[i] = true;
        }
        else
        {
            preset[i] = false;
        }
    }

    if (set() && t2.set())
    {
        direction correspondence[3]{0, 0, 0};
        short signd[3];

        for (direction i=0; i<3; i++)
        {
            if (preset[i])
            {
                signd[i] = 0;
                continue;
            }

            scalar mostAligned = -1;
            for (direction j=0; j<3; j++)
            {
                bool set = false;
                for (direction k=0; k<i; k++)
                {
                    if (correspondence[k] == j)
                    {
                        set = true;
                        break;
                    }
                }

                if (!set)
                {
                    scalar a = operator[](i) & t2.operator[](j);
                    scalar maga = mag(a);

                    if (maga > mostAligned)
                    {
                        correspondence[i] = j;
                        mostAligned = maga;
                        signd[i] = sign(a);
                    }
                }
            }

            operator[](i) += signd[i]*t2.operator[](correspondence[i]);
        }
    }
}


void Foam::triad::align(const vector& v)
{
    if (set())
    {
        vector mostAligned
        (
            mag(v & operator[](0)),
            mag(v & operator[](1)),
            mag(v & operator[](2))
        );

        scalar mav;

        if
        (
            mostAligned.x() > mostAligned.y()
         && mostAligned.x() > mostAligned.z()
        )
        {
            mav = mostAligned.x();
            mostAligned = operator[](0);
        }
        else if (mostAligned.y() > mostAligned.z())
        {
            mav = mostAligned.y();
            mostAligned = operator[](1);
        }
        else
        {
            mav = mostAligned.z();
            mostAligned = operator[](2);
        }

        if (mav < 0.99)
        {
            tensor R(rotationTensor(mostAligned, v));

            operator[](0) = transform(R, operator[](0));
            operator[](1) = transform(R, operator[](1));
            operator[](2) = transform(R, operator[](2));
        }
    }
}


Foam::triad Foam::triad::sortxyz() const
{
    if (!this->set())
    {
        return *this;
    }

    triad t;

    if
    (
        mag(operator[](0).x()) > mag(operator[](1).x())
     && mag(operator[](0).x()) > mag(operator[](2).x())
    )
    {
        t[0] = operator[](0);

        if (mag(operator[](1).y()) > mag(operator[](2).y()))
        {
            t[1] = operator[](1);
            t[2] = operator[](2);
        }
        else
        {
            t[1] = operator[](2);
            t[2] = operator[](1);
        }
    }
    else if
    (
        mag(operator[](1).x()) > mag(operator[](2).x())
    )
    {
        t[0] = operator[](1);

        if (mag(operator[](0).y()) > mag(operator[](2).y()))
        {
            t[1] = operator[](0);
            t[2] = operator[](2);
        }
        else
        {
            t[1] = operator[](2);
            t[2] = operator[](0);
        }
    }
    else
    {
        t[0] = operator[](2);

        if (mag(operator[](0).y()) > mag(operator[](1).y()))
        {
            t[1] = operator[](0);
            t[2] = operator[](1);
        }
        else
        {
            t[1] = operator[](1);
            t[2] = operator[](0);
        }
    }

    if (t[0].x() < 0) t[0] *= -1;
    if (t[1].y() < 0) t[1] *= -1;
    if (t[2].z() < 0) t[2] *= -1;

    return t;
}



Foam::triad::operator Foam::quaternion() const
{
    tensor R;

    R.xx() = x().x();
    R.xy() = y().x();
    R.xz() = z().x();

    R.yx() = x().y();
    R.yy() = y().y();
    R.yz() = z().y();

    R.zx() = x().z();
    R.zy() = y().z();
    R.zz() = z().z();

    return quaternion(R);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::triad::operator=(const tensor& t)
{
    x() = t.x();
    y() = t.y();
    z() = t.z();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::diff(const triad& A, const triad& B)
{
    triad tmpA = A.sortxyz();
    triad tmpB = B.sortxyz();

    scalar sumDifference = 0;

    for (direction dir = 0; dir < 3; dir++)
    {
        if (!tmpA.set(dir) || !tmpB.set(dir))
        {
            continue;
        }

        scalar cosPhi =
            (tmpA[dir] & tmpB[dir])
           /(mag(tmpA[dir])*mag(tmpA[dir]) + small);

        cosPhi = min(max(cosPhi, -1), 1);

        sumDifference += mag(cosPhi - 1);
    }

    return (sumDifference/3);
}


// ************************************************************************* //

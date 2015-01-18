/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "indexedVertex.H"
#include "point.H"
#include "Istream.H"
#include "Ostream.H"
#include "OStringStream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    CGAL::Point_3<baseK>& p
)
{
//    string data(is);
//
//    std::istringstream stdIs;
//
//    CGAL::set_ascii_mode(stdIs);
//
//    stdIs.str(data);
//
//    CGAL::Gmpz xNumer, xDenom;
//    CGAL::Gmpz yNumer, yDenom;
//    CGAL::Gmpz zNumer, zDenom;
//
//    stdIs >> xNumer >> xDenom >> yNumer >> yDenom >> zNumer >> zDenom;
//
//    CGAL::Gmpq x(xNumer, xDenom);
//    CGAL::Gmpq y(yNumer, yDenom);
//    CGAL::Gmpq z(zNumer, zDenom);
//
//    p = CGAL::Point_3<baseK>
//    (
//        CGAL::to_double(x),
//        CGAL::to_double(y),
//        CGAL::to_double(z)
//    );

    Foam::point pt;

    is >> pt.x() >> pt.y() >> pt.z();

    p = CGAL::Point_3<baseK>
    (
        pt.x(),
        pt.y(),
        pt.z()
    );

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CGAL::Point_3<baseK>& p
)
{
//    CGAL::Gmpq x(CGAL::to_double(p.x()));
//    CGAL::Gmpq y(CGAL::to_double(p.y()));
//    CGAL::Gmpq z(CGAL::to_double(p.z()));
//
//    std::ostringstream stdOs;
//
//    CGAL::set_ascii_mode(stdOs);
//
//    stdOs<< x.numerator() << ' ' << x.denominator() << ' '
//         << y.numerator() << ' ' << y.denominator() << ' '
//         << z.numerator() << ' ' << z.denominator();
//
//    os << stdOs.str();

    os  << CGAL::to_double(p.x()) << ' '
        << CGAL::to_double(p.y()) << ' '
        << CGAL::to_double(p.z());

    return os;
}


template<class Gt, class Vb>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CGAL::indexedVertex<Gt, Vb>& p
)
{
    os  << p.point() << ' '
        << p.index() << ' '
        << static_cast<int>(p.type()) << ' '
        << p.procIndex() << ' '
        << p.alignment() << ' '
        << p.targetCellSize() << ' '
        << static_cast<int>(p.fixed());

    return os;
}


template<class Gt, class Vb>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    CGAL::indexedVertex<Gt, Vb>& p
)
{
    is  >> p.point()
        >> p.index();

    int type;
    is  >> type;
    p.type() = static_cast<Foam::indexedVertexEnum::vertexType>(type);

    is  >> p.procIndex()
        >> p.alignment()
        >> p.targetCellSize();

    int fixed;
    is  >> fixed;
    p.fixed() = static_cast<bool>(fixed);

    return is;
}


template<class Gt, class Vb>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<CGAL::indexedVertex<Gt, Vb> >& p
)
{
    const CGAL::indexedVertex<Gt, Vb>& iv = p.t_;

    const Foam::point pt
    (
        CGAL::to_double(iv.point().x()),
        CGAL::to_double(iv.point().y()),
        CGAL::to_double(iv.point().z())
    );

    string fixed
    (
        iv.vertexFixed_
      ? string(" fixed, ")
      : string(" free, ")
    );

    string referred
    (
        Pstream::myProcNo() == iv.processor_
      ? string(" (local)")
      : string(" (from " + name(iv.processor_) + ")")
    );

    os  << iv.index_ << " "
        << CGAL::indexedVertex<Gt, Vb>::vertexTypeNames_[iv.type_]
        << " at:" << pt
        << " size:" << iv.targetCellSize_
        << " alignment:" << iv.alignment_
        << fixed
        << referred.c_str()
        << endl;

    return os;
}


// ************************************************************************* //

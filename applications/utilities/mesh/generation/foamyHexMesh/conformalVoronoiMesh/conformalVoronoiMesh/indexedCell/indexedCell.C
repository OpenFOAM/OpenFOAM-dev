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

#include "indexedCell.H"

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

template<class Gt, class Cb>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<CGAL::indexedCell<Gt, Cb> >& p
)
{
    const CGAL::indexedCell<Gt, Cb>& iv = p.t_;

    os  << "Cell: ";

    if (iv.index_ == CGAL::indexedCell<Gt, Cb>::ctFar)
    {
        os  << "far";
    }
    else if (iv.index_ >= 0)
    {
        os  << iv.index_;
    }
    else if (iv.index_ == CGAL::indexedCell<Gt, Cb>::ctInternal)
    {
        os  << "internal";
    }
    else if (iv.index_ == CGAL::indexedCell<Gt, Cb>::ctSurface)
    {
        os  << "surface";
    }
    else if (iv.index_ == CGAL::indexedCell<Gt, Cb>::ctFeatureEdge)
    {
        os  << "featureEdge";
    }
    else if (iv.index_ == CGAL::indexedCell<Gt, Cb>::ctFeaturePoint)
    {
        os  << "featurePoint";
    }
    else
    {
        os  << "unassigned";
    }

    if (iv.parallelDualVertex())
    {
        os  << " (processor)";
    }
    else
    {
        os  << " (local)";
    }

    os  << " filterCount: " << iv.filterCount_ << nl;
    os  << "    " << iv.vertex(0)->info();
    os  << "    " << iv.vertex(1)->info();
    os  << "    " << iv.vertex(2)->info();
    os  << "    " << iv.vertex(3)->info();

    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Read construct cellPrimitiveModel from Istream.
    Write cellPrimitiveModel to Ostream

\*---------------------------------------------------------------------------*/

#include "cellModel.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cellModel::cellModel(Istream& is)
{
    dictionaryEntry entry(dictionary::null, is);
    name_ = entry.keyword();
    entry.lookup("index") >> index_;
    entry.lookup("numberOfPoints") >> nPoints_;
    entry.lookup("faces") >> faces_;
    entry.lookup("edges") >> edges_;
}


Ostream& operator<<(Ostream& os, const cellModel& c)
{
    os  << "name" << tab << c.name_ << tab
        << "index" << tab << c.index_ << tab
        << "numberOfPoints" << tab << c.nPoints_ << tab
        << "faces" << tab << c.faces_ << tab
        << "edges" << tab << c.edges_ << endl;

    return os;
}


#if defined (__GNUC__)
template<>
#endif
Ostream& operator<<(Ostream& os, const InfoProxy<cellModel>& ip)
{
    const cellModel& cm = ip.t_;

    os  << "name = " << cm.name() << ", "
        << "index = " << cm.index() << ", "
        << "number of points = " << cm.nPoints() << ", "
        << "number of faces = " << cm.nFaces() << ", "
        << "number of edges = " << cm.nEdges()
        << endl;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::face::face()
{}


inline Foam::face::face(label s)
:
    labelList(s, -1)
{}


inline Foam::face::face(const labelUList& lst)
:
    labelList(lst)
{}


inline Foam::face::face(const labelList& lst)
:
    labelList(lst)
{}


inline Foam::face::face(labelList&& lst)
:
    labelList(move(lst))
{}


inline Foam::face::face(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::pointField Foam::face::points(const pointField& meshPoints) const
{
    // There are as many points as there labels for them
    pointField p(size());

    // For each point in list, set it to the point in 'pnts' addressed
    // by 'labs'
    forAll(p, i)
    {
        p[i] = meshPoints[operator[](i)];
    }

    // Return list
    return p;
}


inline Foam::scalar Foam::face::mag(const pointField& points) const
{
    return ::Foam::mag(area(points));
}


inline Foam::label Foam::face::nEdges() const
{
    // for a closed polygon a number of edges is the same as number of points
    return size();
}


inline Foam::edge Foam::face::faceEdge(const label n) const
{
    return edge(operator[](n), operator[](fcIndex(n)));
}


inline Foam::label Foam::face::nextLabel(const label i) const
{
    return operator[](fcIndex(i));
}


inline Foam::label Foam::face::prevLabel(const label i) const
{
    return operator[](rcIndex(i));
}


inline Foam::label Foam::face::nTriangles() const
{
    return size() - 2;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::face::operator=(labelList&& l)
{
    labelList::operator=(move(l));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const face& a, const face& b)
{
    return face::compare(a,b) != 0;
}


inline bool Foam::operator!=(const face& a, const face& b)
{
    return face::compare(a,b) == 0;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

inline Foam::Istream& Foam::operator>>(Istream& is, face& f)
{
    is >> static_cast<labelList&>(f);

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, face&)");

    return is;
}


// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "blockDescriptor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockDescriptor::blockDescriptor
(
    const cellShape& bshape,
    const pointField& blockPointField,
    const curvedEdgeList& edges,
    const Vector<label>& meshDensity,
    const UList<scalar>& expand,
    const word& zoneName
)
:
    blockPointField_(blockPointField),
    curvedEdges_(edges),
    blockShape_(bshape),
    meshDensity_(meshDensity),
    edgePoints_(12),
    edgeWeights_(12),
    expand_(expand),
    zoneName_(zoneName)
{
    if (expand_.size() != 12)
    {
        FatalErrorIn
        (
            "blockDescriptor::blockDescriptor"
            "(const cellShape&, const pointField& blockPointField, "
            "const curvedEdgeList&, const Vector<label>& meshDensity, "
            "const scalarList& expand, const word& zoneName)"
        )   << "Unknown definition of expansion ratios"
            << exit(FatalError);
    }

    // create a list of edges
    makeBlockEdges();
}


Foam::blockDescriptor::blockDescriptor
(
    const pointField& blockPointField,
    const curvedEdgeList& edges,
    Istream& is
)
:
    blockPointField_(blockPointField),
    curvedEdges_(edges),
    blockShape_(is),
    meshDensity_(),
    edgePoints_(12),
    edgeWeights_(12),
    expand_(12, 1.0),
    zoneName_()
{
    // Examine next token
    token t(is);

    // Optional zone name
    if (t.isWord())
    {
        zoneName_ = t.wordToken();

        // Examine next token
        is >> t;
    }
    is.putBack(t);

    if (t.isPunctuation())
    {
        // new-style: read a list of 3 values
        if (t.pToken() == token::BEGIN_LIST)
        {
            is >> meshDensity_;
        }
        else
        {
            FatalIOErrorIn
            (
                "blockDescriptor::blockDescriptor"
                "(const pointField&, const curvedEdgeList&, Istream&)",
                is
            )   << "incorrect token while reading n, expected '(', found "
                << t.info()
                << exit(FatalIOError);
        }
    }
    else
    {
        // old-style: read three labels
        is  >> meshDensity_.x()
            >> meshDensity_.y()
            >> meshDensity_.z();
    }

    is >> t;
    if (!t.isWord())
    {
        is.putBack(t);
    }

    scalarList expRatios(is);

    if (expRatios.size() == 1)
    {
        // identical in x/y/z-directions
        expand_ = expRatios[0];
    }
    else if (expRatios.size() == 3)
    {
        // x-direction
        expand_[0]  = expRatios[0];
        expand_[1]  = expRatios[0];
        expand_[2]  = expRatios[0];
        expand_[3]  = expRatios[0];

        // y-direction
        expand_[4]  = expRatios[1];
        expand_[5]  = expRatios[1];
        expand_[6]  = expRatios[1];
        expand_[7]  = expRatios[1];

        // z-direction
        expand_[8]  = expRatios[2];
        expand_[9]  = expRatios[2];
        expand_[10] = expRatios[2];
        expand_[11] = expRatios[2];
    }
    else if (expRatios.size() == 12)
    {
        expand_ = expRatios;
    }
    else
    {
        FatalErrorIn
        (
            "blockDescriptor::blockDescriptor"
            "(const pointField&, const curvedEdgeList&, Istream&)"
        )   << "Unknown definition of expansion ratios: " << expRatios
            << exit(FatalError);
    }

    // create a list of edges
    makeBlockEdges();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockDescriptor::~blockDescriptor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::blockDescriptor::blockPointField() const
{
    return blockPointField_;
}


const Foam::cellShape& Foam::blockDescriptor::blockShape() const
{
    return blockShape_;
}


const Foam::List< Foam::List< Foam::point > >&
Foam::blockDescriptor::blockEdgePoints() const
{
    return edgePoints_;
}


const Foam::scalarListList& Foam::blockDescriptor::blockEdgeWeights() const
{
    return edgeWeights_;
}


const Foam::Vector<Foam::label>& Foam::blockDescriptor::meshDensity() const
{
    return meshDensity_;
}


const Foam::word& Foam::blockDescriptor::zoneName() const
{
    return zoneName_;
}


Foam::label Foam::blockDescriptor::nPoints() const
{
    return
    (
        (meshDensity_.x() + 1)
      * (meshDensity_.y() + 1)
      * (meshDensity_.z() + 1)
    );
}


Foam::label Foam::blockDescriptor::nCells() const
{
    return
    (
        meshDensity_.x()
      * meshDensity_.y()
      * meshDensity_.z()
    );
}


const Foam::point& Foam::blockDescriptor::blockPoint(const label i) const
{
    return blockPointField_[blockShape_[i]];
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const blockDescriptor& bd)
{
    const cellShape& bshape = bd.blockShape();
    const labelList& blockLabels = bshape;

    os  << bshape.model().name() << " (";

    forAll(blockLabels, labelI)
    {
        if (labelI)
        {
            os  << ' ';
        }
        os  << blockLabels[labelI];
    }
    os  << ')';

    if (bd.zoneName().size())
    {
        os  << ' ' << bd.zoneName();
    }

    os  << ' '  << bd.meshDensity()
        << " simpleGrading (";


    const scalarList& expand = bd.expand_;

    // can we use a compact notation?
    if
    (
        // x-direction
        (
            expand[0] == expand[1]
         && expand[0] == expand[2]
         && expand[0] == expand[3]
        )
     && // y-direction
        (
            expand[4] == expand[5]
         && expand[4] == expand[6]
         && expand[4] == expand[7]
        )
     && // z-direction
        (
            expand[8] == expand[9]
         && expand[8] == expand[10]
         && expand[8] == expand[11]
        )
    )
    {
        os  << expand[0] << ' ' << expand[4] << ' ' << expand[8];
    }
    else
    {
        forAll(expand, edgeI)
        {
            if (edgeI)
            {
                os  << ' ';
            }
            os  << expand[edgeI];
        }
    }


    os  << ")";

    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "cellEdgeAddressing.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::cellEdgeAddressing::cellEdgeAddressing
(
    const cellEdgeAddressingList& cAddrs,
    const label celli
)
:
    cAddrs_(cAddrs),
    celli_(celli)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::CompactListList<Foam::label>&
Foam::cellEdgeAddressingData::cfiAndFeiToCei() const
{
    return cfiAndFeiToCei_;
}


inline const Foam::List<Foam::Pair<Foam::labelPair>>&
Foam::cellEdgeAddressingData::ceiToCfiAndFei() const
{
    return ceiToCfiAndFei_;
}


inline const Foam::boolList& Foam::cellEdgeAddressingData::cOwns() const
{
    return cOwns_;
}


inline const Foam::CompactListList<Foam::label>&
Foam::cellEdgeAddressing::cfiAndFeiToCei() const
{
    return cAddrs_.data(celli_).cfiAndFeiToCei();
}


inline const Foam::List<Foam::Pair<Foam::labelPair>>&
Foam::cellEdgeAddressing::ceiToCfiAndFei() const
{
    return cAddrs_.data(celli_).ceiToCfiAndFei();
}


inline const Foam::boolList& Foam::cellEdgeAddressing::cOwns() const
{
    return cAddrs_.data(celli_).cOwns();
}


inline const Foam::cellEdgeAddressingData&
Foam::cellEdgeAddressingList::data(const label celli) const
{
    if (!list_.set(celli))
    {
        const cell& c = mesh().cells()[celli];
        const faceList& fs = mesh().faces();
        const labelList& fOwners = mesh().faceOwner();

        list_.set
        (
            celli,
            new cellEdgeAddressingData(c, fs, fOwners[c[0]] == celli)
        );
    }

    return list_[celli];
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline const Foam::cellEdgeAddressing
Foam::cellEdgeAddressingList::operator[](const label celli) const
{
    return cellEdgeAddressing(*this, celli);
}


// ************************************************************************* //

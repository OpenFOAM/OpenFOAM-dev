/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "multiRegionPrefixer.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::multiRegionPrefixer::prefixes() const
{
    return prefixSingleRegion_ || regionNames_.size() > 1;
}


Foam::string::size_type Foam::multiRegionPrefixer::prefixWidth() const
{
    string::size_type n = 0;

    if (prefixes())
    {
        forAll(regionNames_, regionj)
        {
            n = max(n, regionNames_[regionj].size() + 1);
        }
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiRegionPrefixer::regionPrefixer::regionPrefixer
(
    const multiRegionPrefixer& mrp,
    const label regioni
)
:
    mrp_(mrp),
    regioni_(regioni)
{
    if (mrp_.prefixes() && regioni_ != -1)
    {
        const string::size_type dn =
            mrp_.prefixWidth()
          - mrp_.regionNames_[regioni].size();

        Sout.prefix() =
            mrp_.regionNames_[regioni] + string(dn, ' ');
    }
}


Foam::multiRegionPrefixer::regionPrefixer::regionPrefixer
(
    regionPrefixer&& rp
)
:
    regionPrefixer(rp.mrp_, rp.regioni_)
{
    rp.regioni_ = -1;
}


Foam::multiRegionPrefixer::multiRegionPrefixer
(
    const bool prefixSingleRegion,
    const wordList& regionNames
)
:
    prefixSingleRegion_(prefixSingleRegion),
    regionNames_(regionNames)
{
    if (prefixes())
    {
        Sout.prefix() = string(prefixWidth(), ' ');
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::multiRegionPrefixer::regionPrefixer::~regionPrefixer()
{
    if (mrp_.prefixes() && regioni_ != -1)
    {
        Sout.prefix() = string(mrp_.prefixWidth(), ' ');
    }
}


Foam::multiRegionPrefixer::~multiRegionPrefixer()
{
    if (prefixes())
    {
        Sout.prefix() = string::null;
    }
}


// ************************************************************************* //

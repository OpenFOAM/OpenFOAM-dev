/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "MultiRegionRefs.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Region>
bool Foam::MultiRegionRefs<Region>::prefixes() const
{
    return regions_.size() > 1;
}


template<class Region>
Foam::string::size_type Foam::MultiRegionRefs<Region>::prefixWidth() const
{
    string::size_type n = 0;

    if (prefixes())
    {
        forAll(regions_, regioni)
        {
            n = max(n, regionName(regions_[regioni]).size() + 1);
        }
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Region>
template<class NonConstRegion>
Foam::RegionRef<Region>::RegionRef
(
    const MultiRegionRefs<NonConstRegion>& mrr,
    const label regioni,
    Region& region
)
:
    mrr_(*reinterpret_cast<const MultiRegionRefs<Region>*>(&mrr)),
    regioni_(regioni),
    region_(region)
{
    if (mrr_.prefixes() && regioni_ != -1)
    {
        const string::size_type dn =
            mrr_.prefixWidth()
          - regionName(mrr_.regions_[regioni]).size();

        Sout.prefix() =
            regionName(mrr_.regions_[regioni]) + string(dn, ' ');
    }
}


template<class Region>
Foam::RegionRef<Region>::RegionRef(RegionRef&& rp)
:
    RegionRef(rp.mrr_, rp.regioni_, rp.region_)
{
    rp.regioni_ = -1;
}


template<class Region>
Foam::MultiRegionRefs<Region>::MultiRegionRefs(UPtrList<Region>& regions)
:
    regions_(regions),
    previousPrefix_(Sout.prefix()),
    indices_()
{
    forAll(regions, regioni)
    {
        indices_.insert(regionName(regions_[regioni]), regioni);
    }

    if (prefixes())
    {
        Sout.prefix() = string(prefixWidth(), ' ');
    }
    else
    {
        Sout.prefix() = word::null;
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class Region>
Foam::RegionRef<Region>::~RegionRef()
{
    if (mrr_.prefixes() && regioni_ != -1)
    {
        Sout.prefix() = string(mrr_.prefixWidth(), ' ');
    }
}


template<class Region>
Foam::MultiRegionRefs<Region>::~MultiRegionRefs()
{
    Sout.prefix() = previousPrefix_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Region>
Foam::label Foam::MultiRegionRefs<Region>::size() const
{
    return regions_.size();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Region>
Foam::RegionRef<Region>::operator Region&() const
{
    return region_;
}


template<class Region>
Region& Foam::RegionRef<Region>::operator()() const
{
    return region_;
}


template<class Region>
Foam::RegionRef<const Region> Foam::MultiRegionRefs<Region>::operator[]
(
    const label regioni
) const
{
    return RegionRef<const Region>(*this, regioni, regions_[regioni]);
}


template<class Region>
Foam::RegionRef<Region> Foam::MultiRegionRefs<Region>::operator[]
(
    const label regioni
)
{
    return RegionRef<Region>(*this, regioni, regions_[regioni]);
}


template<class Region>
Foam::RegionRef<const Region> Foam::MultiRegionRefs<Region>::operator[]
(
    const word& regionName
) const
{
    return operator[](indices_[regionName]);
}


template<class Region>
Foam::RegionRef<Region> Foam::MultiRegionRefs<Region>::operator[]
(
    const word& regionName
)
{
    return operator[](indices_[regionName]);
}


// ************************************************************************* //

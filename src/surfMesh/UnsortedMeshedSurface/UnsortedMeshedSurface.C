/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canReadType
(
    const word& ext,
    const bool verbose
)
{
   return fileFormats::surfaceFormatsCore::checkSupport
   (
       readTypes() | ParentType::readTypes(),
       ext,
       verbose,
       "reading"
   );
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes(),
        ext,
        verbose,
        "writing"
    );
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canRead
(
    const fileName& name,
    const bool verbose
)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return canReadType(ext, verbose);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::write
(
    const fileName& name,
    const UnsortedMeshedSurface<Face>& surf
)
{
    if (debug)
    {
        InfoInFunction << "Writing to " << name << endl;
    }

    const word ext = name.ext();

    typename writefileExtensionMemberFunctionTable::iterator mfIter =
        writefileExtensionMemberFunctionTablePtr_->find(ext);

    if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
    {
        // no direct writer, delegate to proxy if possible
        wordHashSet supported = ProxyType::writeTypes();

        if (supported.found(ext))
        {
            MeshedSurfaceProxy<Face>(surf).write(name);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown file extension " << ext << nl << nl
                << "Valid types are :" << endl
                << (supported | writeTypes())
                << exit(FatalError);
        }
    }
    else
    {
        mfIter()(name, surf);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface()
:
    ParentType()
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const Xfer<List<label>>& zoneIds,
    const Xfer<surfZoneIdentifierList>& zoneTofc
)
:
    ParentType(pointLst, faceLst),
    zoneIds_(zoneIds),
    zoneToc_(zoneTofc)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const labelUList& zoneSizes,
    const UList<word>& zoneNames
)
:
    ParentType(pointLst, faceLst)
{
    if (zoneSizes.size())
    {
        if (zoneNames.size())
        {
            setZones(zoneSizes, zoneNames);
        }
        else
        {
            setZones(zoneSizes);
        }
    }
    else
    {
        setOneZone();
    }
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
:
    ParentType
    (
        xferCopy(surf.points()),
        xferCopy(surf.faces())
    ),
    zoneIds_(surf.zoneIds()),
    zoneToc_(surf.zoneToc())
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const MeshedSurface<Face>& surf
)
:
    ParentType
    (
        xferCopy(surf.points()),
        xferCopy(surf.faces())
    )
{
    setZones(surf.surfZones());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<UnsortedMeshedSurface<Face>>& surf
)
:
    ParentType()
{
    transfer(surf());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<MeshedSurface<Face>>& surf
)
:
    ParentType()
{
    transfer(surf());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const fileName& name,
    const word& ext
)
:
    ParentType()
{
    read(name, ext);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface(const fileName& name)
:
    ParentType()
{
    read(name);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Time& t,
    const word& surfName
)
:
    ParentType()
{
    MeshedSurface<Face> surf(t, surfName);
    transfer(surf);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::UnsortedMeshedSurface<Face>::~UnsortedMeshedSurface()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setOneZone()
{
    zoneIds_.setSize(size());
    zoneIds_ = 0;

    word zoneName;
    if (zoneToc_.size())
    {
        zoneName = zoneToc_[0].name();
    }
    if (zoneName.empty())
    {
        zoneName = "zone0";
    }

    // set single default zone
    zoneToc_.setSize(1);
    zoneToc_[0] = surfZoneIdentifier(zoneName, 0);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const surfZoneList& zoneLst
)
{
    zoneIds_.setSize(size());
    zoneToc_.setSize(zoneLst.size());

    forAll(zoneToc_, zoneI)
    {
        const surfZone& zone = zoneLst[zoneI];
        zoneToc_[zoneI] = zone;

        // assign sub-zone Ids
        SubList<label> subZone(zoneIds_, zone.size(), zone.start());
        subZone = zoneI;
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const labelUList& sizes,
    const UList<word>& names
)
{
    zoneIds_.setSize(size());
    zoneToc_.setSize(sizes.size());

    label start = 0;
    forAll(zoneToc_, zoneI)
    {
        zoneToc_[zoneI] = surfZoneIdentifier(names[zoneI], zoneI);

        // assign sub-zone Ids
        SubList<label> subZone(zoneIds_, sizes[zoneI], start);
        subZone = zoneI;

        start += sizes[zoneI];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const labelUList& sizes
)
{
    zoneIds_.setSize(size());
    zoneToc_.setSize(sizes.size());

    label start = 0;
    forAll(zoneToc_, zoneI)
    {
        zoneToc_[zoneI] = surfZoneIdentifier
        (
            word("zone") + ::Foam::name(zoneI),
            zoneI
        );

        // assign sub-zone Ids
        SubList<label> subZone(zoneIds_, sizes[zoneI], start);
        subZone = zoneI;

        start += sizes[zoneI];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::remapFaces
(
    const labelUList& faceMap
)
{
    // re-assign the zone Ids
    if (notNull(faceMap) && faceMap.size())
    {
        if (zoneToc_.empty())
        {
            setOneZone();
        }
        else if (zoneToc_.size() == 1)
        {
            // optimized for single-zone case
            zoneIds_ = 0;
        }
        else
        {
            List<label> newZones(faceMap.size());

            forAll(faceMap, facei)
            {
                newZones[facei] = zoneIds_[faceMap[facei]];
            }
            zoneIds_.transfer(newZones);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setSize(const label s)
{
    this->storedFaces().setSize(s);
    // if zones extend: set with last zoneId
    zoneIds_.setSize(s, zoneToc_.size() - 1);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::clear()
{
    ParentType::clear();
    zoneIds_.clear();
    zoneToc_.clear();
}


template<class Face>
Foam::surfZoneList Foam::UnsortedMeshedSurface<Face>::sortedZones
(
    labelList& faceMap
) const
{
    // supply some zone names
    Map<word> zoneNames;
    forAll(zoneToc_, zoneI)
    {
        zoneNames.insert(zoneI, zoneToc_[zoneI].name());
    }

    // std::sort() really seems to mix up the order.
    // and std::stable_sort() might take too long / too much memory

    // Assuming that we have relatively fewer zones compared to the
    // number of items, just do it ourselves

    // step 1: get zone sizes and store (origId => zoneI)
    Map<label> lookup;
    forAll(zoneIds_, facei)
    {
        const label origId = zoneIds_[facei];

        Map<label>::iterator fnd = lookup.find(origId);
        if (fnd != lookup.end())
        {
            fnd()++;
        }
        else
        {
            lookup.insert(origId, 1);
        }
    }

    // step 2: assign start/size (and name) to the newZones
    // re-use the lookup to map (zoneId => zoneI)
    surfZoneList zoneLst(lookup.size());
    label start = 0;
    label zoneI = 0;
    forAllIter(Map<label>, lookup, iter)
    {
        label origId = iter.key();

        word name;
        Map<word>::const_iterator fnd = zoneNames.find(origId);
        if (fnd != zoneNames.end())
        {
            name = fnd();
        }
        else
        {
            name = word("zone") + ::Foam::name(zoneI);
        }

        zoneLst[zoneI] = surfZone
        (
            name,
            0,           // initialize with zero size
            start,
            zoneI
        );

        // increment the start for the next zone
        // and save the (zoneId => zoneI) mapping
        start += iter();
        iter() = zoneI++;
    }


    // step 3: build the re-ordering
    faceMap.setSize(zoneIds_.size());

    forAll(zoneIds_, facei)
    {
        label zoneI = lookup[zoneIds_[facei]];
        faceMap[facei] = zoneLst[zoneI].start() + zoneLst[zoneI].size()++;
    }

    // with reordered faces registered in faceMap
    return zoneLst;
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>
Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField&  locPoints = this->localPoints();
    const List<Face>&  locFaces  = this->localFaces();

    // Fill pointMap, faceMap
    PatchTools::subsetMap(*this, include, pointMap, faceMap);

    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList  oldToNew(locPoints.size());
    forAll(pointMap, pointi)
    {
        newPoints[pointi] = locPoints[pointMap[pointi]];
        oldToNew[pointMap[pointi]] = pointi;
    }

    // Renumber face node labels and compact
    List<Face>  newFaces(faceMap.size());
    List<label> newZones(faceMap.size());

    forAll(faceMap, facei)
    {
        const label origFacei = faceMap[facei];
        newFaces[facei] = Face(locFaces[origFacei]);

        // Renumber labels for face
        Face& f = newFaces[facei];
        forAll(f, fp)
        {
            f[fp] = oldToNew[f[fp]];
        }

        newZones[facei] = zoneIds_[origFacei];
    }
    oldToNew.clear();

    // construct a sub-surface
    return UnsortedMeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newZones),
        xferCopy(zoneToc_)
    );
}


template<class Face>
Foam::UnsortedMeshedSurface<Face> Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include
) const
{
    labelList pointMap, faceMap;
    return subsetMesh(include, pointMap, faceMap);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::reset
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const Xfer<List<label>>& zoneIds
)
{
    ParentType::reset
    (
        pointLst,
        faceLst,
        Xfer<surfZoneList>()
    );

    if (notNull(zoneIds))
    {
        zoneIds_.transfer(zoneIds());
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::reset
(
    const Xfer<List<point>>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const Xfer<List<label>>& zoneIds
)
{
    ParentType::reset
    (
        pointLst,
        faceLst,
        Xfer<surfZoneList>()
    );

    if (notNull(zoneIds))
    {
        zoneIds_.transfer(zoneIds());
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    UnsortedMeshedSurface<Face>& surf
)
{
    ParentType::reset
    (
        xferMove(surf.storedPoints()),
        xferMove(surf.storedFaces()),
        Xfer<surfZoneList>()
    );

    zoneIds_.transfer(surf.zoneIds_);
    zoneToc_.transfer(surf.zoneToc_);

    surf.clear();
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    ParentType::reset
    (
        xferMove(surf.storedPoints()),
        xferMove(surf.storedFaces()),
        Xfer<surfZoneList>()
    );

    setZones(surf.surfZones());
    surf.clear();
}


template<class Face>
Foam::Xfer<Foam::UnsortedMeshedSurface<Face>>
Foam::UnsortedMeshedSurface<Face>::xfer()
{
    return xferMove(*this);
}


// Read from file, determine format from extension
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();
        return read(unzipName, unzipName.ext());
    }
    else
    {
        return read(name, ext);
    }
}


// Read from file in given format
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read
(
    const fileName& name,
    const word& ext
)
{
    clear();

    // read via use selector mechanism
    transfer(New(name, ext)());
    return true;
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::write
(
    const Time& t,
    const word& surfName
) const
{
    MeshedSurfaceProxy<Face>(*this).write(t, surfName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::operator=
(
    const UnsortedMeshedSurface<Face>& surf
)
{
    clear();

    this->storedPoints() = surf.points();
    this->storedFaces()  = surf.faces();
    zoneIds_ = surf.zoneIds_;
    zoneToc_ = surf.zoneToc_;
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::operator
Foam::MeshedSurfaceProxy<Face>() const
{
    labelList faceMap;
    List<surfZone> zoneLst = this->sortedZones(faceMap);

    return MeshedSurfaceProxy<Face>
    (
        this->points(),
        this->faces(),
        zoneLst,
        faceMap
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UnsortedMeshedSurfaceNew.C"

// ************************************************************************* //

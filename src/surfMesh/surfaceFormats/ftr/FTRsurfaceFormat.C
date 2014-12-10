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

#include "FTRsurfaceFormat.H"
#include "Keyed.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::FTRsurfaceFormat<Face>::FTRsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::FTRsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::FTRsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    List<ftrPatch> ftrPatches(is);

    // points read directly
    is >> this->storedPoints();

    // triFaces read with attached keys
    List< Keyed<triFace> > facesRead(is);

    List<Face>  faceLst(facesRead.size());
    List<label> zoneIds(facesRead.size());

    // disentangle faces/keys - already triangulated
    forAll(facesRead, faceI)
    {
        // unfortunately cannot transfer to save memory
        faceLst[faceI] = facesRead[faceI];
        zoneIds[faceI] = facesRead[faceI].key();
    }

    this->storedFaces().transfer(faceLst);
    this->storedZoneIds().transfer(zoneIds);
    facesRead.clear();

    // change ftrPatch into surfZoneIdentifier
    List<surfZoneIdentifier> newZones(ftrPatches.size());
    forAll(newZones, zoneI)
    {
        newZones[zoneI] = surfZoneIdentifier
        (
            ftrPatches[zoneI].name(),
            zoneI
        );
    }

    this->storedZoneToc().transfer(newZones);
    return true;
}


// ************************************************************************* //

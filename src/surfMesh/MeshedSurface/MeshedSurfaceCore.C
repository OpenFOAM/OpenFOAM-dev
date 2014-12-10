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

#include "MeshedSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

    // specialization: from face -> triFace
    template<>
    void Foam::MeshedSurface<triFace>::transcribe(MeshedSurface<face>& surf)
    {
        // first triangulate
        surf.triangulate();
        this->storedPoints().transfer(surf.storedPoints());
        this->storedZones().transfer(surf.storedZones());

        // transcribe from face -> triFace
        List<face>&    origFaces = surf.storedFaces();
        List<triFace>  newFaces(origFaces.size());
        forAll(origFaces, faceI)
        {
            newFaces[faceI] = triFace
            (
                static_cast<const labelUList&>(origFaces[faceI])
            );
        }
        surf.clear();

        this->storedFaces().transfer(newFaces);
    }


    // specialization: from face -> face
    template<>
    void Foam::MeshedSurface<face>::transcribe(MeshedSurface<face>& surf)
    {
        this->transfer(surf);
    }


}  // end of namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

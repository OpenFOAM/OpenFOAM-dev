/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "triSurface.H"
#include "STLtriangle.H"
#include "primitivePatch.H"
#include "HashTable.H"
#include "hashSignedLabel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::writeSTLASCII(const bool writeSorted, Ostream& os) const
{
    labelList faceMap;

    surfacePatchList patches(calcPatches(faceMap));

    if (writeSorted)
    {
        label faceIndex = 0;
        forAll(patches, patchi)
        {
            // Print all faces belonging to this region
            const surfacePatch& patch = patches[patchi];

            os  << "solid " << patch.name() << endl;

            for
            (
                label patchFacei = 0;
                patchFacei < patch.size();
                patchFacei++
            )
            {
                const label facei = faceMap[faceIndex++];

                const vector& n = faceNormals()[facei];

                os  << "  facet normal "
                    << n.x() << ' ' << n.y() << ' ' << n.z() << nl
                    << "    outer loop" << endl;

                const labelledTri& f = (*this)[facei];
                const point& pa = points()[f[0]];
                const point& pb = points()[f[1]];
                const point& pc = points()[f[2]];

                os  << "       vertex "
                    << pa.x() << ' ' << pa.y() << ' ' << pa.z() << nl
                    << "       vertex "
                    << pb.x() << ' ' << pb.y() << ' ' << pb.z() << nl
                    << "       vertex "
                    << pc.x() << ' ' << pc.y() << ' ' << pc.z() << nl
                    << "    endloop" << nl
                    << "  endfacet" << endl;
            }

            os  << "endsolid " << patch.name() << endl;
        }
    }
    else
    {
        // Get patch (=compact region) per face
        labelList patchIDs(size());
        forAll(patches, patchi)
        {
            label facei = patches[patchi].start();

            forAll(patches[patchi], i)
            {
                patchIDs[faceMap[facei++]] = patchi;
            }
        }

        label currentPatchi = -1;

        forAll(*this, facei)
        {
            if (currentPatchi != patchIDs[facei])
            {
                if (currentPatchi != -1)
                {
                    // Have already valid patch. Close it.
                    os  << "endsolid " << patches[currentPatchi].name()
                        << nl;
                }
                currentPatchi = patchIDs[facei];
                os  << "solid " << patches[currentPatchi].name() << nl;
            }

            const vector& n = faceNormals()[facei];

            os  << "  facet normal "
                << n.x() << ' ' << n.y() << ' ' << n.z() << nl
                << "    outer loop" << endl;

            const labelledTri& f = (*this)[facei];
            const point& pa = points()[f[0]];
            const point& pb = points()[f[1]];
            const point& pc = points()[f[2]];

            os  << "       vertex "
                << pa.x() << ' ' << pa.y() << ' ' << pa.z() << nl
                << "       vertex "
                << pb.x() << ' ' << pb.y() << ' ' << pb.z() << nl
                << "       vertex "
                << pc.x() << ' ' << pc.y() << ' ' << pc.z() << nl
                << "    endloop" << nl
                << "  endfacet" << endl;
        }

        if (currentPatchi != -1)
        {
            os  << "endsolid " << patches[currentPatchi].name()
                << nl;
        }
    }
}


void Foam::triSurface::writeSTLBINARY(std::ostream& os) const
{
    // Write the STL header
    string header("Foam binary STL");
    header.resize(STLheaderSize);
    os.write(header.c_str(), STLheaderSize);

    label nTris = size();
    os.write(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    const vectorField& normals = faceNormals();

    forAll(*this, facei)
    {
        const labelledTri& f = (*this)[facei];

        // Convert vector into STL single precision
        STLpoint n(normals[facei]);
        STLpoint pa(points()[f[0]]);
        STLpoint pb(points()[f[1]]);
        STLpoint pc(points()[f[2]]);

        STLtriangle stlTri(n, pa, pb, pc, f.region());

        stlTri.write(os);
    }
}


// ************************************************************************* //

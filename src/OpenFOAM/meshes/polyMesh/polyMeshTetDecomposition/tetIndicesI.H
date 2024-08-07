/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "tetIndices.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::tetIndices::cell() const
{
    return celli_;
}


inline Foam::label& Foam::tetIndices::cell()
{
    return celli_;
}


inline Foam::label Foam::tetIndices::face() const
{
    return facei_;
}


inline Foam::label& Foam::tetIndices::face()
{
    return facei_;
}


inline Foam::label Foam::tetIndices::tetPt() const
{
    return tetPti_;
}


inline Foam::label& Foam::tetIndices::tetPt()
{
    return tetPti_;
}


inline Foam::triFace Foam::tetIndices::faceTriIs(const polyMesh& mesh) const
{
    const Foam::face& f = mesh.faces()[face()];

    label faceBasePtI = mesh.tetBasePtIs()[face()];

    if (faceBasePtI < 0)
    {
        static labelHashSet badFaces;
        static label badTimeIndex = -1;

        if (badTimeIndex != mesh.time().timeIndex())
        {
            badFaces.clear();
            badTimeIndex = mesh.time().timeIndex();
        }

        if (!badFaces[face()])
        {
            WarningInFunction
                << "No base point for face " << face() << ", " << f
                << ", produces a valid tet decomposition." << endl;

            badFaces.set(face());
        }

        faceBasePtI = 0;
    }

    label facePtI = (tetPt() + faceBasePtI) % f.size();
    label faceOtherPtI = f.fcIndex(facePtI);

    if (mesh.faceOwner()[face()] != cell())
    {
        Swap(facePtI, faceOtherPtI);
    }

    return triFace(f[faceBasePtI], f[facePtI], f[faceOtherPtI]);
}


inline Foam::tetPointRef Foam::tetIndices::tet
(
    const polyMesh& mesh,
    const pointField& meshPoints,
    const pointField& cellCentres
) const
{
    const triFace tri = faceTriIs(mesh);

    return tetPointRef
    (
        cellCentres[cell()],
        meshPoints[tri[0]],
        meshPoints[tri[1]],
        meshPoints[tri[2]]
    );
}


inline Foam::tetPointRef Foam::tetIndices::tet(const polyMesh& mesh) const
{
    return tet(mesh, mesh.points(), mesh.cellCentres());
}


inline Foam::triPointRef Foam::tetIndices::faceTri
(
    const polyMesh& mesh,
    const pointField& meshPoints
) const
{
    const triFace tri = faceTriIs(mesh);

    return triPointRef
    (
        meshPoints[tri[0]],
        meshPoints[tri[1]],
        meshPoints[tri[2]]
    );
}


inline Foam::triPointRef Foam::tetIndices::faceTri(const polyMesh& mesh) const
{
    return faceTri(mesh, mesh.points());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::tetIndices::operator==(const Foam::tetIndices& rhs) const
{
    return
        cell() == rhs.cell()
     && face() == rhs.face()
     && tetPt() == rhs.tetPt();
}


inline bool Foam::tetIndices::operator!=(const Foam::tetIndices& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //

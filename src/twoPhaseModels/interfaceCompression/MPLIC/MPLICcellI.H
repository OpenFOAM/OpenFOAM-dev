/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "MPLICcell.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::MPLICcell::calcAddressing
(
    const MPLICcellStorage& cellInfo
)
{
    localEdgeFaces_.setSize(cellInfo.cellEdges().size());
    localFaceEdges_.setSize(cellInfo.size());

    // Create edge map, setting edge faces to zero
    Map<label> edgeMap;
    forAll(cellInfo.cellEdges(), edgei)
    {
        edgeMap.set(cellInfo.cellEdges()[edgei], edgei);
        localEdgeFaces_[edgei].clear();
    }

    // Set size of each face edges to zero
    forAll(localFaceEdges_, i)
    {
        localFaceEdges_[i].clear();
    }

    forAll(cellInfo.cellFaces(), facei)
    {
        const label faceN = cellInfo.cellFaces()[facei];
        const labelList& faceEdges = cellInfo.faceEdges()[faceN];
        forAll(faceEdges, edgei)
        {
            const label edg = faceEdges[edgei];
            const label localEdgeIndex = edgeMap[edg];
            localEdgeFaces_[localEdgeIndex].append(facei);
            localFaceEdges_[facei].append(localEdgeIndex);
        }
    }
    addressingCalculated_ = true;
}


inline bool Foam::MPLICcell::cutStatusCalcSf()
{
    bool cutOrientationDiffers = false;
    const point pAvg = sum(cutPoints_)/cutPoints_.size();
    for (label i=0; i<cutPoints_.size(); i+=2)
    {
        const vector area = (cutPoints_[i] - pAvg)^(cutPoints_[i+1] - pAvg);
        cutSf_ += area;
        if
        (
            sign(area.x()*cutSf_.x()) == -1
         || sign(area.y()*cutSf_.y()) == -1
         || sign(area.z()*cutSf_.z()) == -1
        )
        {
            cutOrientationDiffers = true;
        }
    }
    cutSf_ *= 0.5;

    return cutOrientationDiffers;
}


inline Foam::vector Foam::MPLICcell::calcCutSf() const
{
    vector cutArea = Zero;
    const point pAvg = sum(cutPoints_)/cutPoints_.size();
    for (label i=0; i<cutPoints_.size(); i+=2)
    {
        cutArea += (cutPoints_[i] - pAvg)^(cutPoints_[i+1] - pAvg);
    }
    cutArea *= 0.5;

    return cutArea;
}


inline Foam::vector Foam::MPLICcell::calcCutCf(const vector& cutSf) const
{
    const vector sumAHat = normalised(cutSf);
    const point pAvg = sum(cutPoints_)/cutPoints_.size();

    scalar sumAn = 0;
    vector sumAnc = Zero;

    for (label i=0; i < cutPoints_.size(); i+=2)
    {
        const vector a = (cutPoints_[i] - pAvg) ^ (cutPoints_[i+1] - pAvg);
        const vector c = cutPoints_[i] + cutPoints_[i+1] + pAvg;
        const scalar an = a & sumAHat;
        sumAn += an;
        sumAnc += an*c;
    }

    if (sumAn > vSmall)
    {
        return (1.0/3.0)*sumAnc/sumAn;
    }
    else
    {
        return pAvg;
    }
}


inline void Foam::MPLICcell::appendSfCf
(
    const vector& Sf,
    const vector& Cf,
    const scalar magSf,
    const bool own
)
{
    if (magSf > 0)
    {
        if (own)
        {
            subFaceAreas_.append(Sf);
        }
        else
        {
            subFaceAreas_.append(-Sf);
        }
        subFaceCentres_.append(Cf);
    }
}


inline void Foam::MPLICcell::phiU
(
    const pointField& points,
    const faceList& faces,
    const labelList& cFaces,
    const vectorField& pointsU
)
{
    const label nFaces = cFaces.size();

    // Set size and initialise alphaPhiU
    alphaPhiU_.setSize(nFaces);
    alphaPhiU_ = 0;

    // Set size and initialise phiU
    phiU_.setSize(nFaces);

    // Reconstruct fluxes
    forAll(cFaces, facei)
    {
        phiU_[facei] =
            MPLICface().alphaPhiU(pointsU, points, faces[cFaces[facei]]);
    }
}


inline void Foam::MPLICcell::resetFaceFields(const label nFaces)
{
    alphaf_.setSize(nFaces);
    alphaf_ = 0;
    if (unweighted_)
    {
        subFaceMagSf_.setSize(nFaces);
        subFaceMagSf_ = 0;
    }
    else
    {
        alphaPhiU_.setSize(nFaces);
        alphaPhiU_ = 0;
    }
}


inline void Foam::MPLICcell::calcAlphaf
(
    const UIndirectList<scalar>& magSfs
)
{
    const label nFaces = magSfs.size();
    alphaf_.setSize(nFaces);
    forAll(alphaf_, facei)
    {
        if (magSfs[facei] > vSmall)
        {
            alphaf_[facei] = subFaceMagSf_[facei]/magSfs[facei];

            // Bound between 0 and 1 (it is always > 0)
            alphaf_[facei] = (alphaf_[facei] > 1) ? 1 : alphaf_[facei];
        }
        else
        {
            alphaf_[facei] = 0;
        }
    }
}


inline void Foam::MPLICcell::calcAlphaUf()
{
    const label nFaces = alphaPhiU_.size();
    alphaf_.setSize(nFaces);
    forAll(alphaf_, facei)
    {
        if (mag(phiU_[facei]) > vSmall)
        {
            alphaf_[facei] = alphaPhiU_[facei]/phiU_[facei];

            // Bound between 0 and 1
            alphaf_[facei] = (alphaf_[facei] > 1) ? 1 : alphaf_[facei];
            alphaf_[facei] = (alphaf_[facei] < 0) ? 0 : alphaf_[facei];
        }
        else
        {
            alphaf_[facei] = 0;
        }
    }
}


inline void Foam::MPLICcell::clearOneCut()
{
    cutPoints_.clear();
    cutEdges_.clear();
    subFaceAreas_.clear();
    subFaceCentres_.clear();
}


inline void Foam::MPLICcell::clear()
{
    clearOneCut();
    alphaPhiU_.clear();
    subFaceMagSf_.clear();
    cutAlpha_ = -1;
    cutNormal_ = Zero;
    cutSf_ = Zero;
    subCellVolume_ = 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::DynamicList<Foam::scalar>& Foam::MPLICcell::alphaf() const
{
    return alphaf_;
}


inline const Foam::vector& Foam::MPLICcell::cutNormal() const
{
    return cutNormal_;
}


inline Foam::scalar Foam::MPLICcell::subCellVolume() const
{
    return subCellVolume_;
}


inline Foam::scalar Foam::MPLICcell::cutAlpha() const
{
    return cutAlpha_;
}


// ************************************************************************* //

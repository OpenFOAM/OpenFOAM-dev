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

#include "regionModel1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::labelListList&
Foam::regionModels::regionModel1D::boundaryFaceFaces() const
{
    return boundaryFaceFaces_;
}


inline const Foam::labelListList&
Foam::regionModels::regionModel1D::boundaryFaceCells() const
{
    return boundaryFaceCells_;
}


inline const Foam::labelList&
Foam::regionModels::regionModel1D::boundaryFaceOppositeFace() const
{
    return boundaryFaceOppositeFace_;
}


inline const Foam::surfaceScalarField&
Foam::regionModels::regionModel1D::nMagSf() const
{
    if (!nMagSfPtr_.valid())
    {
        FatalErrorIn
        (
            "inline const Foam::surfaceScalarField&"
            "Foam::regionModel1Ds::regionModel1D::nMagSf() const"
        )<< "Face normal areas not available" << abort(FatalError);
    }

    return nMagSfPtr_();
}


inline Foam::label Foam::regionModels::regionModel1D::nLayers() const
{
    return nLayers_;
}


// ************************************************************************* //

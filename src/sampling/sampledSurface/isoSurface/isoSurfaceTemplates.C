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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurface::interpolate
(
    const Field<Type>& cellCoords,
    const Field<Type>& pointCoords
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(pointToVerts_.size()));
    Field<Type>& fld = tfld.ref();

    forAll(pointToVerts_, i)
    {
        scalar s0;
        Type p0;
        {
            label v0 = pointToVerts_[i][0];
            if (v0 < mesh_.nPoints())
            {
                s0 = pVals_[v0];
                p0 = pointCoords[v0];
            }
            else
            {
                label celli = v0-mesh_.nPoints();
                s0 = cVals_[celli];
                p0 = cellCoords[celli];
            }
        }

        scalar s1;
        Type p1;
        {
            label v1 = pointToVerts_[i][1];
            if (v1 < mesh_.nPoints())
            {
                s1 = pVals_[v1];
                p1 = pointCoords[v1];
            }
            else
            {
                label celli = v1-mesh_.nPoints();
                s1 = cVals_[celli];
                p1 = cellCoords[celli];
            }
        }

        scalar d = s1-s0;
        if (mag(d) > VSMALL)
        {
            scalar s = (iso_-s0)/d;
            fld[i] = s*p1+(1.0-s)*p0;
        }
        else
        {
            fld[i] = 0.5*(p0+p1);
        }
    }

    return tfld;
}


// ************************************************************************* //

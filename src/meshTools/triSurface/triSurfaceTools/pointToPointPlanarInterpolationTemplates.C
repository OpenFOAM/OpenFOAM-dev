/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "pointToPointPlanarInterpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::pointToPointPlanarInterpolation::interpolate
(
    const Field<Type>& sourceFld
) const
{
    if (nPoints_ != sourceFld.size())
    {
        FatalErrorInFunction
            << "Number of source points = " << nPoints_
            << " number of values = " << sourceFld.size()
            << exit(FatalError);
    }

    tmp<Field<Type>> tfld(new Field<Type>(nearestVertex_.size()));
    Field<Type>& fld = tfld.ref();

    forAll(fld, i)
    {
        const FixedList<label, 3>& verts = nearestVertex_[i];
        const FixedList<scalar, 3>& w = nearestVertexWeight_[i];

        if (verts[2] == -1)
        {
            if (verts[1] == -1)
            {
                // Use vertex0 only
                fld[i] = sourceFld[verts[0]];
            }
            else
            {
                // Use vertex 0,1
                fld[i] =
                    w[0]*sourceFld[verts[0]]
                  + w[1]*sourceFld[verts[1]];
            }
        }
        else
        {
            fld[i] =
                w[0]*sourceFld[verts[0]]
              + w[1]*sourceFld[verts[1]]
              + w[2]*sourceFld[verts[2]];
        }
    }
    return tfld;
}


// ************************************************************************* //

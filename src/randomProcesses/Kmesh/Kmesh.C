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

#include "Kmesh.H"
#include "polyMesh.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::label Foam::Kmesh::index
(
    const label i,
    const label j,
    const label k,
    const labelList& nn
)
{
    return (k + j*nn[2] + i*nn[1]*nn[2]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Kmesh::Kmesh(const fvMesh& mesh)
:
    vectorField(mesh.V().size()),
    nn_(vector::dim)
{
    boundBox box = mesh.bounds();
    l_ = box.span();

    vector cornerCellCentre = ::Foam::max(mesh.C().primitiveField());
    vector cellL = 2*(box.max() - cornerCellCentre);

    vector rdeltaByL;
    label nTot = 1;

    forAll(nn_, i)
    {
        nn_[i] = label(l_[i]/cellL[i] + 0.5);
        nTot *= nn_[i];

        if (nn_[i] > 1)
        {
            l_[i] -= cellL[i];
        }

        rdeltaByL[i] = nn_[i]/(l_[i]*l_[i]);
    }

    if (nTot != mesh.nCells())
    {
        FatalErrorInFunction
            << "calculated number of cells is incorrect"
            << abort(FatalError);
    }

    for (label i=0; i<nn_[0]; i++)
    {
        scalar k1 = (i - nn_[0]/2)*constant::mathematical::twoPi/l_[0];

        for (label j=0; j<nn_[1]; j++)
        {
            scalar k2 = (j - nn_[1]/2)*constant::mathematical::twoPi/l_[1];

            for (label k=0; k<nn_[2]; k++)
            {
                scalar k3 = (k - nn_[2]/2)*constant::mathematical::twoPi/l_[2];

                (*this)[index(i, j, k, nn_)] = vector(k1, k2, k3);
            }
        }
    }

    kmax_ = mag
    (
        Foam::max
        (
            cmptMag((*this)[index(nn_[0]-1, nn_[1]-1, nn_[2]-1, nn_)]),
            cmptMag((*this)[index(0, 0, 0, nn_)])
        )
    );
}


// ************************************************************************* //

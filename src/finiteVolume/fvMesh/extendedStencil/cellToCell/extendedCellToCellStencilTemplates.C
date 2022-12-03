/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "extendedCellToCellStencil.H"
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class WeightType>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<WeightType, Type>::type>
> Foam::extendedCellToCellStencil::weightedSum
(
    const distributionMap& map,
    const labelListList& stencil,
    const VolField<Type>& fld,
    const List<List<WeightType>>& stencilWeights
)
{
    typedef typename outerProduct<WeightType, Type>::type WeightedType;

    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type>> stencilFld;
    extendedCellToFaceStencil::collectData(map, stencil, fld, stencilFld);

    tmp<VolField<WeightedType>> twf
    (
        new VolField<WeightedType>
        (
            IOobject
            (
                fld.name(),
                mesh.time().name(),
                mesh
            ),
            mesh,
            dimensioned<WeightedType>
            (
                fld.name(),
                fld.dimensions(),
                Zero
            )
        )
    );
    VolField<WeightedType>& wf = twf();

    forAll(wf, celli)
    {
        const List<Type>& stField = stencilFld[celli];
        const List<WeightType>& stWeight = stencilWeights[celli];

        forAll(stField, i)
        {
            wf[celli] += stWeight[i]*stField[i];
        }
    }

    // Boundaries values?

    return twf;
}


// ************************************************************************* //

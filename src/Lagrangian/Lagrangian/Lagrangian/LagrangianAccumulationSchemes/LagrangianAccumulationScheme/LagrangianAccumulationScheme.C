/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianAccumulationScheme.H"
#include "LagrangianSubMesh.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Lagrangian::accumulationScheme<Type>>
Foam::Lagrangian::accumulationScheme<Type>::New
(
    const LagrangianMesh& mesh,
    Istream& is
)
{
    if (is.eof())
    {
        FatalIOErrorInFunction(is)
            << "Accumulation scheme not specified" << endl << endl
            << "Valid accumulation schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(is);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(is)
            << "Unknown accumulation scheme " << schemeName << nl << nl
            << "Valid accumulation schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Lagrangian::accumulationScheme<Type>::~accumulationScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
template<class CellMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, CellMesh>>
Foam::Lagrangian::accumulationScheme<Type>::accumulate
(
    const DimensionedField<Type, LagrangianMesh, PrimitiveField>& lPsi
)
{
    typedef typename DimensionedField<Type, CellMesh>::Mesh resultMeshType;

    tmp<DimensionedField<Type, CellMesh>> tResult =
        DimensionedField<Type, CellMesh>::New
        (
            lPsi.name(),
            refCast<const resultMeshType>(lPsi.mesh().mesh()),
            dimensioned<Type>(lPsi.dimensions(), pTraits<Type>::zero)
        );

    accumulate
    (
        lPsi.mesh().subAll().sub(lPsi),
        tResult.ref().primitiveFieldRef()
    );

    return tResult;
}


template<class Type>
template<class CellMesh, template<class> class PrimitiveField>
void Foam::Lagrangian::accumulationScheme<Type>::accumulate
(
    const LagrangianSubField<Type, PrimitiveField>& lPsi,
    DimensionedField<Type, CellMesh>& vPsi
)
{
    vPsi.dimensions() += lPsi.dimensions();

    accumulate(toSubField(lPsi)(), vPsi.primitiveFieldRef());
}


// ************************************************************************* //

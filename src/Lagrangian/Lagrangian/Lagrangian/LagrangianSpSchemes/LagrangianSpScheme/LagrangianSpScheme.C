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

#include "LagrangianSpScheme.H"
#include "LagrangianSubFields.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class Type, class SpType>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::Lagrangian::SpScheme<Type, SpType>::inner
(
    const LagrangianSubScalarField& Sp,
    const LagrangianSubSubField<Type>& psi
)
{
    return Sp*psi;
}


template<class Type, class SpType>
Foam::tmp
<
    Foam::LagrangianSubField
    <
        typename Foam::innerProduct<Foam::tensor, Type>::type
    >
>
Foam::Lagrangian::SpScheme<Type, SpType>::inner
(
    const LagrangianSubTensorField& Sp,
    const LagrangianSubSubField<Type>& psi
)
{
    return Sp & psi;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type, class SpType>
Foam::tmp<Foam::Lagrangian::SpScheme<Type, SpType>>
Foam::Lagrangian::SpScheme<Type, SpType>::New
(
    const LagrangianMesh& mesh,
    Istream& is
)
{
    if (is.eof())
    {
        FatalIOErrorInFunction(is)
            << "Sp scheme not specified" << endl << endl
            << "Valid Sp schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(is);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(is)
            << "Unknown Sp scheme " << schemeName << nl << nl
            << "Valid Sp schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class SpType>
Foam::Lagrangian::SpScheme<Type, SpType>::~SpScheme()
{}


// ************************************************************************* //

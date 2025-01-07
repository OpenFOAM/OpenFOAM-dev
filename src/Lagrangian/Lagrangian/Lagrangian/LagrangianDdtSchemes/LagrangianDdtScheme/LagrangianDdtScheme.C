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

#include "LagrangianDdtScheme.H"
#include "LagrangianSubFields.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Lagrangian::ddtScheme<Type>>
Foam::Lagrangian::ddtScheme<Type>::New
(
    const LagrangianMesh& mesh,
    Istream& is
)
{
    if (is.eof())
    {
        FatalIOErrorInFunction(is)
            << "Ddt scheme not specified" << endl << endl
            << "Valid ddt schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(is);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(is)
            << "Unknown ddt scheme " << schemeName << nl << nl
            << "Valid ddt schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Lagrangian::ddtScheme<Type>::~ddtScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtScheme<Type>::Lagrangianmddt
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubSubField<Type>& psi
)
{
    tmp<LagrangianEqn<Type>> tEqn(new LagrangianEqn<Type>(deltaT, psi));

    tEqn.ref().deltaTSp += 1;
    tEqn.ref().deltaTSu -= psi.oldTime();

    return tEqn;
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtScheme<Type>::Lagrangianmddt
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    tmp<LagrangianEqn<Type>> tEqn(new LagrangianEqn<Type>(deltaT, psi));

    tEqn.ref().deltaTSp += m;
    tEqn.ref().deltaTSu -= m.oldTime()*psi.oldTime();

    return tEqn;
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtScheme<Type>::Lagrangianmddt0
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubSubField<Type>& psi
)
{
    tmp<LagrangianEqn<Type>> tEqn(new LagrangianEqn<Type>(deltaT, psi));

    tEqn.ref().deltaTSu += psi - psi.oldTime();

    return tEqn;
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtScheme<Type>::Lagrangianmddt0
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    tmp<LagrangianEqn<Type>> tEqn(new LagrangianEqn<Type>(deltaT, psi));

    tEqn.ref().deltaTSu += m*psi - m.oldTime()*psi.oldTime();

    return tEqn;
}


// ************************************************************************* //

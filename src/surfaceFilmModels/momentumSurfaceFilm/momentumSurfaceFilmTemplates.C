/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "momentumSurfaceFilm.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Type> Foam::momentumSurfaceFilm::constrainFilmField
(
    const tmp<Type>& tfield,
    const typename Type::cmptType& value
)
{
    tmp<Type> tresult(tfield);
    Type& result = tresult.ref();

    typename Type::Boundary& fieldBf = result.boundaryFieldRef();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        fieldBf[patchi] = value;

        DebugInFunction
            << "Constraining " << tfield().name()
            << " boundary " << tfield().boundaryField()[patchi].patch().name()
            << " to " << value << endl;
    }

    forAll(passivePatchIDs(), i)
    {
        const label patchi = passivePatchIDs()[i];
        fieldBf[patchi] = value;
        DebugInFunction
            << "Constraining " << tfield().name()
            << " boundary " << tfield().boundaryField()[patchi].patch().name()
            << " to " << value << endl;
    }

    return tresult;
}


// ************************************************************************* //

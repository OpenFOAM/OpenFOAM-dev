/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

InNamespace
    Foam::fvi

Description
    Calculate the gradient of the given field returning a volInternalField

SourceFiles
    fviGrad.C

\*---------------------------------------------------------------------------*/


#ifndef fviGrad_H
#define fviGrad_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvi functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvi
{
    template<class Type>
    tmp<VolInternalField<typename outerProduct<vector, Type>::type>>
    grad
    (
        const SurfaceField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<typename outerProduct<vector, Type>::type>>
    grad
    (
        const tmp<SurfaceField<Type>>&
    );

    template<class Type>
    tmp<VolInternalField<typename outerProduct<vector, Type>::type>>
    grad
    (
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<VolInternalField<typename outerProduct<vector, Type>::type>>
    grad
    (
        const tmp<VolField<Type>>&,
        const word& name
    );

    template<class Type>
    tmp<VolInternalField<typename outerProduct<vector, Type>::type>>
    grad
    (
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<typename outerProduct<vector, Type>::type>>
    grad
    (
        const tmp<VolField<Type>>&
    );

    template<class Expression, class = EnableIfExpression<Expression>>
    inline auto grad(const Expression& e)
    {
        return grad(eval(e));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fviGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

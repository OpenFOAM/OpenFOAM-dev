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
    Surface integrate surfaceField  returning a volInternalField.
    Surface sum a surfaceField  returning a volInternalField.

SourceFiles
    fviSurfaceIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef fviSurfaceIntegrate_H
#define fviSurfaceIntegrate_H

#include "primitiveFieldsFwd.H"
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
    void surfaceIntegrate
    (
        Field<Type>&,
        const SurfaceField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>>
    surfaceIntegrate(const SurfaceField<Type>&);

    template<class Type>
    tmp<VolInternalField<Type>>
    surfaceIntegrate(const tmp<SurfaceField<Type>>&);

    template<class Expression, class = EnableIfExpression<Expression>>
    auto surfaceIntegrate(const Expression&);

    template<class Type>
    tmp<VolInternalField<Type>> surfaceSum(const SurfaceField<Type>&);

    template<class Type>
    tmp<VolInternalField<Type>> surfaceSum(const tmp<SurfaceField<Type>>&);

    template<class Expression, class = EnableIfExpression<Expression>>
    auto surfaceSum(const Expression&);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fviSurfaceIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

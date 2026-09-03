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
    Calculate the laplacian of the given field returning a volInternalField

SourceFiles
    fviLaplacian.C

\*---------------------------------------------------------------------------*/


#ifndef fviLaplacian_H
#define fviLaplacian_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvi functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvi
{
    template<class Type>
    tmp<VolInternalField<Type>> laplacian
    (
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Expression, class = EnableIfExpression<Expression>>
    inline auto laplacian
    (
        const Expression& e,
        const word& name
    )
    {
        return laplacian(eval(e), name);
    }

    template<class Type>
    tmp<VolInternalField<Type>> laplacian
    (
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<VolField<Type>>&
    );

    template<class Expression, class = EnableIfExpression<Expression>>
    inline auto laplacian
    (
        const Expression& e
    )
    {
        return laplacian(eval(e));
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template
    <
        class GType,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const dimensioned<GType>& d,
        const Expression& e,
        const word& name
    )
    {
        return laplacian(d, eval(e), name);
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const tmp<VolField<Type>>&
    );

    template
    <
        class GType,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const dimensioned<GType>& dt,
        const Expression& e
    )
    {
        return laplacian(dt, eval(e));
    }


    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const VolField<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template
    <
        class GType,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const VolField<GType>& d,
        const Expression& e,
        const word& name
    )
    {
        return laplacian(d, eval(e), name);
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const VolField<GType>&,
        const tmp<VolField<Type>>&
    );

    template
    <
        class GType,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const VolField<GType>& d,
        const Expression& e
    )
    {
        return laplacian(d, eval(e));
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const tmp<VolField<Type>>&
    );


    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template
    <
        class Type,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline tmp<VolInternalField<Type>> laplacian
    (
        const Expression& e,
        const VolField<Type>& vf,
        const word& name
    )
    {
        return laplacian(eval(e), vf, name);
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template
    <
        class GType,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const SurfaceField<GType>& d,
        const Expression& e,
        const word& name
    )
    {
        return laplacian(d, eval(e), name);
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template
    <
        class DExpression,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const DExpression& d,
        const Expression& e,
        const word& name
    )
    {
        return laplacian(eval(d), eval(e), name);
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&
    );

    template
    <
        class Type,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline tmp<VolInternalField<Type>> laplacian
    (
        const Expression& e,
        const VolField<Type>& vf
    )
    {
        return laplacian(eval(e), vf);
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const tmp<VolField<Type>>&
    );

    template
    <
        class GType,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const SurfaceField<GType>& d,
        const Expression& e
    )
    {
        return laplacian(d, eval(e));
    }

    template<class Type, class GType>
    tmp<VolInternalField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const tmp<VolField<Type>>&
    );

    template
    <
        class DExpression,
        class Expression,
        class = EnableIfExpression<Expression>
    >
    inline auto laplacian
    (
        const DExpression& d,
        const Expression& e
    )
    {
        return laplacian(eval(d), eval(e));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fviLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

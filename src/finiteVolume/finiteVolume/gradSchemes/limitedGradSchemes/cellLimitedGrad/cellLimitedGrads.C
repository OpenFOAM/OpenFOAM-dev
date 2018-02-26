/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "cellLimitedGrad.H"
#include "minmodGradientLimiter.H"
#include "VenkatakrishnanGradientLimiter.H"
#include "cubicGradientLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeNamedFvLimitedGradTypeScheme(SS, Type, Limiter, Name)              \
    typedef Foam::fv::SS<Foam::Type, Foam::fv::gradientLimiters::Limiter>      \
        SS##Type##Limiter##_;                                                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Type##Limiter##_,                                                  \
        Name,                                                                  \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            gradScheme<Type>::addIstreamConstructorToTable                     \
            <                                                                  \
                SS<Type, gradientLimiters::Limiter>                            \
            > add##SS##Type##Limiter##IstreamConstructorToTable_;              \
        }                                                                      \
    }

#define makeFvLimitedGradTypeScheme(SS, Type, Limiter)                         \
    makeNamedFvLimitedGradTypeScheme(SS##Grad, Type, Limiter, #SS"<"#Limiter">")

#define makeFvLimitedGradScheme(SS, Limiter)                                   \
                                                                               \
    makeFvLimitedGradTypeScheme(SS, scalar, Limiter)                           \
    makeFvLimitedGradTypeScheme(SS, vector, Limiter)


// Default limiter in minmod specified without the limiter name
// for backward compatibility
makeNamedFvLimitedGradTypeScheme(cellLimitedGrad, scalar, minmod, "cellLimited")
makeNamedFvLimitedGradTypeScheme(cellLimitedGrad, vector, minmod, "cellLimited")

makeFvLimitedGradScheme(cellLimited, Venkatakrishnan)
makeFvLimitedGradScheme(cellLimited, cubic)

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

Class
    Foam::interpolations::pointMVC

Description
    Uses volPointInterpolation to create values on the points. Then does a Mean
    Value Coordinates (MVC) interpolation.

SourceFiles
    pointMVCI.H
    pointMVC.C

\*---------------------------------------------------------------------------*/

#ifndef pointMVC_H
#define pointMVC_H

#include "volPointInterpolation_interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolations
{

/*---------------------------------------------------------------------------*\
                         Class pointMVC Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointMVC
:
    public fieldInterpolation<Type, pointMVC<Type>>,
    public volPointInterpolation<Type>
{
public:

    //- Runtime type information
    TypeName("pointMVC");


    // Constructors

        //- Construct from components
        pointMVC(const VolField<Type>& psi);

        //- Copy construct
        pointMVC(const pointMVC<Type>& i);

        //- Clone
        virtual autoPtr<interpolation<Type>> clone() const
        {
            return autoPtr<interpolation<Type>>
            (
                new pointMVC<Type>(*this)
            );
        }


    // Member Functions

        //- Inherit interpolate from interpolation
        using interpolation<Type>::interpolate;

        //- Interpolate field to the given point in the given cell
        inline Type interpolate
        (
            const vector& position,
            const label celli,
            const label facei = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointMVCI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointMVC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

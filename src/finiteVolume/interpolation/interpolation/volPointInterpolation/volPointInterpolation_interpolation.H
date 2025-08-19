/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::interpolations::volPointInterpolation

Description
    Base class for interpolations that require a vol-point interpolated field

\*---------------------------------------------------------------------------*/

#ifndef volPointInterpolation_interpolation_H
#define volPointInterpolation_interpolation_H

#include "interpolation.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolations
{

/*---------------------------------------------------------------------------*\
                   Class volPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class volPointInterpolation
{
protected:

    // Protected data

        //- Interpolated volfield
        const PointField<Type> psip_;


public:

    //- Runtime type information
    TypeName("volPointInterpolation");


    // Constructors

        //- Construct from components
        volPointInterpolation(const VolField<Type>& psi);

        //- Copy construct
        volPointInterpolation(const volPointInterpolation<Type>& i);

        //- Construct from components
        volPointInterpolation
        (
            const VolField<Type>& psi,
            tmp<PointField<Type>> psip
        );


    // Destructor
    virtual ~volPointInterpolation();


    // Member Functions

        //- Access the point field
        inline const PointField<Type>& psip() const
        {
            return psip_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volPointInterpolation_interpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

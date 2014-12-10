/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    pointLinear

Description
    Face-point interpolation scheme class derived from linear and
    returns linear weighting factors but also applies an explicit correction.

    Uses volPointInterpolation to obtain the field values at the face-points.

SourceFiles
    pointLinear.C

\*---------------------------------------------------------------------------*/

#ifndef pointLinear_H
#define pointLinear_H

#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pointLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointLinear
:
    public linear<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointLinear(const pointLinear&);

        //- Disallow default bitwise assignment
        void operator=(const pointLinear&);


public:

    //- Runtime type information
    TypeName("pointLinear");


    // Constructors

        //- Construct from mesh
        pointLinear(const fvMesh& mesh)
        :
            linear<Type>(mesh)
        {}


        //- Construct from mesh and Istream
        pointLinear
        (
            const fvMesh& mesh,
            Istream&
        )
        :
            linear<Type>(mesh)
        {}


        //- Construct from mesh, faceFlux and Istream
        pointLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            linear<Type>(mesh)
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

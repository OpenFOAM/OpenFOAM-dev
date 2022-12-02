/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::deferred

Description
    Deferred correction interpolation scheme derived from upwind
    which returns upwind weighting factors and an explicit correction obtained
    from the specified scheme.

    This ensures that the transport matrix generated is at least diagonally
    equal and avoids the need for relaxation of the matrix (which can affect
    conservation) for transient running.

SourceFiles
    deferred.C

\*---------------------------------------------------------------------------*/

#ifndef deferred_H
#define deferred_H

#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class deferred Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class deferred
:
    public upwind<Type>
{
    // Private member data

        //- Scheme
        tmp<surfaceInterpolationScheme<Type>> tScheme_;


public:

    //- Runtime type information
    TypeName("deferred");


    // Constructors

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        deferred
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, schemeData),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, schemeData)
            )
        {}

        //- Construct from faceFlux and Istream
        deferred
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, faceFlux, schemeData),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New
                (
                    mesh,
                    faceFlux,
                    schemeData
                )
            )
        {}

        //- Disallow default bitwise copy construction
        deferred(const deferred&) = delete;


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>& vf
        ) const
        {
            return SurfaceField<Type>::New
            (
                "deferred::correction(" + vf.name() + ')',
                tScheme_->interpolate(vf)
              - upwind<Type>::interpolate(vf, upwind<Type>::weights())
            );
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const deferred&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

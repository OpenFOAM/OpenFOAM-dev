/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::fv::gaussGrad

Description
    Basic second-order gradient scheme using face-interpolation
    and Gauss' theorem.

SourceFiles
    gaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef gaussGrad_H
#define gaussGrad_H

#include "gradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class gaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussGrad
:
    public fv::gradScheme<Type>
{

protected:

    // Protected Data

        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;


    // Protected Member Functions

        static void calcGrad
        (
            VolInternalField<typename outerProduct<vector, Type>::type>& grad,
            const SurfaceField<Type>& sf
        );


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from mesh
        gaussGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linear<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        gaussGrad(const fvMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(nullptr)
        {
            if (is.eof())
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type>>
                    (
                        new linear<Type>(mesh)
                    );
            }
            else
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type>>
                    (
                        surfaceInterpolationScheme<Type>::New(mesh, is)
                    );
            }
        }

        //- Disallow default bitwise copy construction
        gaussGrad(const gaussGrad&) = delete;


    // Member Functions

        //- Set the grad field to the gradient of the volField vf
        virtual void calcGrad
        (
            VolInternalField<typename outerProduct<vector, Type>::type>& grad,
            const VolField<Type>& vf
        ) const;

        // using fv::gradScheme<Type>::fviGrad;
        using fv::gradScheme<Type>::fvcGrad;

        //- Return the gradient of the given field
        //  calculated using Gauss' theorem on the given surface field
        static tmp<VolInternalField<typename outerProduct<vector, Type>::type>>
        fviGrad
        (
            const SurfaceField<Type>&,
            const word& name
        );

        //- Return the gradient of the given field
        //  calculated using Gauss' theorem on the given surface field
        static tmp<VolField<typename outerProduct<vector, Type>::type>>
        fvcGrad
        (
            const SurfaceField<Type>&,
            const word& name
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gaussGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

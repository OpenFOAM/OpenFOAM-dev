/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::fv::laplacianScheme

Description
    Abstract base class for laplacian schemes.

SourceFiles
    laplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef laplacianScheme_H
#define laplacianScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "linear.H"
#include "correctedSnGrad.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class laplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class laplacianScheme
:
    public tmp<laplacianScheme<Type, GType>>::refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;

        tmp<surfaceInterpolationScheme<GType>> tinterpGammaScheme_;

        tmp<snGradScheme<Type>> tsnGradScheme_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            laplacianScheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        laplacianScheme(const fvMesh& mesh)
        :
            mesh_(mesh),
            tinterpGammaScheme_(new linear<GType>(mesh)),
            tsnGradScheme_(new correctedSnGrad<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        laplacianScheme(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            tinterpGammaScheme_(nullptr),
            tsnGradScheme_(nullptr)
        {
            tinterpGammaScheme_ = tmp<surfaceInterpolationScheme<GType>>
            (
                surfaceInterpolationScheme<GType>::New(mesh, is)
            );

            tsnGradScheme_ = tmp<snGradScheme<Type>>
            (
                snGradScheme<Type>::New(mesh, is)
            );
        }

        //- Construct from mesh, interpolation and snGradScheme schemes
        laplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            mesh_(mesh),
            tinterpGammaScheme_(igs),
            tsnGradScheme_(sngs)
        {}

        //- Disallow default bitwise copy construction
        laplacianScheme(const laplacianScheme&) = delete;


    // Selectors

        //- Return a pointer to a new laplacianScheme created on freestore
        static tmp<laplacianScheme<Type, GType>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~laplacianScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<fvMatrix<Type>> fvmLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        ) = 0;

        virtual tmp<fvMatrix<Type>> fvmLaplacian
        (
            const VolField<GType>&,
            const VolField<Type>&
        );

        virtual tmp<VolField<Type>> fvcLaplacian
        (
            const VolField<Type>&
        ) = 0;

        virtual tmp<VolField<Type>> fvcLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        ) = 0;

        virtual tmp<VolField<Type>> fvcLaplacian
        (
            const VolField<GType>&,
            const VolField<Type>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laplacianScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvLaplacianTypeScheme(SS, GType, Type)                             \
    typedef Foam::fv::SS<Foam::Type, Foam::GType> SS##Type##GType;             \
    defineNamedTemplateTypeNameAndDebug(SS##Type##GType, 0);                   \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            typedef SS<Type, GType> SS##Type##GType;                           \
                                                                               \
            laplacianScheme<Type, GType>::                                     \
                addIstreamConstructorToTable<SS<Type, GType>>                  \
                add##SS##Type##GType##IstreamConstructorToTable_;              \
        }                                                                      \
    }


#define makeFvLaplacianScheme(SS)                                              \
                                                                               \
makeFvLaplacianTypeScheme(SS, scalar, scalar)                                  \
makeFvLaplacianTypeScheme(SS, symmTensor, scalar)                              \
makeFvLaplacianTypeScheme(SS, tensor, scalar)                                  \
makeFvLaplacianTypeScheme(SS, scalar, vector)                                  \
makeFvLaplacianTypeScheme(SS, symmTensor, vector)                              \
makeFvLaplacianTypeScheme(SS, tensor, vector)                                  \
makeFvLaplacianTypeScheme(SS, scalar, sphericalTensor)                         \
makeFvLaplacianTypeScheme(SS, symmTensor, sphericalTensor)                     \
makeFvLaplacianTypeScheme(SS, tensor, sphericalTensor)                         \
makeFvLaplacianTypeScheme(SS, scalar, symmTensor)                              \
makeFvLaplacianTypeScheme(SS, symmTensor, symmTensor)                          \
makeFvLaplacianTypeScheme(SS, tensor, symmTensor)                              \
makeFvLaplacianTypeScheme(SS, scalar, tensor)                                  \
makeFvLaplacianTypeScheme(SS, symmTensor, tensor)                              \
makeFvLaplacianTypeScheme(SS, tensor, tensor)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "laplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

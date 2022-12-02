/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::fv::snGradScheme

Description
    Abstract base class for snGrad schemes.

SourceFiles
    snGradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef snGradScheme_H
#define snGradScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class snGradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class snGradScheme
:
    public tmp<snGradScheme<Type>>::refCount
{
    // Private Data

        //- Hold reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            snGradScheme,
            Mesh,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        snGradScheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        snGradScheme(const snGradScheme&);


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<snGradScheme<Type>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~snGradScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }


        //- Return the snGrad of the given cell field with the given deltaCoeffs
        static tmp<SurfaceField<Type>>
        snGrad
        (
            const VolField<Type>&,
            const tmp<surfaceScalarField>&,
            const word& snGradName = "snGrad"
        );

        //- Return the sndGrad of the given cell field
        static tmp<SurfaceField<Type>>
        sndGrad
        (
            const VolField<Type>&,
            const word& snGradName = "sndGrad"
        );

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>&
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the snGrad
        //  for the given field
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>&) const
        {
            return tmp<SurfaceField<Type>>
            (
                nullptr
            );
        }

        //- Return the snGrad of the given cell field
        //  with explicit correction
        virtual tmp<SurfaceField<Type>>
        snGrad(const VolField<Type>&) const;

        //- Return the snGrad of the given tmp cell field
        //  with explicit correction
        tmp<SurfaceField<Type>>
        snGrad
        (
            const tmp<VolField<Type>>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const snGradScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeSnGradTypeScheme(SS, Type)                                         \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            snGradScheme<Type>::addMeshConstructorToTable<SS<Type>>            \
                add##SS##Type##MeshConstructorToTable_;                        \
        }                                                                      \
    }

#define makeSnGradScheme(SS)                                                   \
                                                                               \
makeSnGradTypeScheme(SS, scalar)                                               \
makeSnGradTypeScheme(SS, vector)                                               \
makeSnGradTypeScheme(SS, sphericalTensor)                                      \
makeSnGradTypeScheme(SS, symmTensor)                                           \
makeSnGradTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "snGradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

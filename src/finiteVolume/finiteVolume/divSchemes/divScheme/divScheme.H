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
    Foam::fv::divScheme

Description
    Abstract base class for div schemes.

SourceFiles
    divScheme.C

\*---------------------------------------------------------------------------*/

#ifndef divScheme_H
#define divScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "linear.H"
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
                           Class divScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class divScheme
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;
        tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;


    // Private Member Functions

        //- Disallow copy construct
        divScheme(const divScheme&);

        //- Disallow default bitwise assignment
        void operator=(const divScheme&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            divScheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        divScheme(const fvMesh& mesh)
        :
            mesh_(mesh),
            tinterpScheme_(new linear<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        divScheme(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            tinterpScheme_(surfaceInterpolationScheme<Type>::New(mesh, is))
        {}


    // Selectors

        //- Return a pointer to a new divScheme created on freestore
        static tmp<divScheme<Type> > New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~divScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp
        <
            GeometricField
            <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
        > fvcDiv
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvDivTypeScheme(SS, Type)                                          \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
divScheme<Type>::addIstreamConstructorToTable<SS<Type> >                       \
    add##SS##Type##IstreamConstructorToTable_;


#define makeFvDivScheme(SS)                                                    \
                                                                               \
makeFvDivTypeScheme(SS, vector)                                                \
makeFvDivTypeScheme(SS, sphericalTensor)                                       \
makeFvDivTypeScheme(SS, symmTensor)                                            \
makeFvDivTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "divScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

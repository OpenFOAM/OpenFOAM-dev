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
    Foam::fv::d2dt2Scheme

Description
    Abstract base class for d2dt2 schemes.

SourceFiles
    d2dt2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef d2dt2Scheme_H
#define d2dt2Scheme_H

#include "tmp.H"
#include "dimensionedType.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
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
                           Class d2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class d2dt2Scheme
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        d2dt2Scheme(const d2dt2Scheme&);

        //- Disallow default bitwise assignment
        void operator=(const d2dt2Scheme&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            d2dt2Scheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        d2dt2Scheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from mesh and Istream
        d2dt2Scheme(const fvMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new d2dt2Scheme created on freestore
        static tmp<d2dt2Scheme<Type> > New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~d2dt2Scheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<fvMatrix<Type> > fvmD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<fvMatrix<Type> > fvmD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<fvMatrix<Type> > fvmD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvD2dt2TypeScheme(SS, Type)                                        \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
d2dt2Scheme<Type>::addIstreamConstructorToTable<SS<Type> >                     \
    add##SS##Type##IstreamConstructorToTable_;


#define makeFvD2dt2Scheme(SS)                                                  \
                                                                               \
makeFvD2dt2TypeScheme(SS, scalar)                                              \
makeFvD2dt2TypeScheme(SS, vector)                                              \
makeFvD2dt2TypeScheme(SS, sphericalTensor)                                     \
makeFvD2dt2TypeScheme(SS, symmTensor)                                          \
makeFvD2dt2TypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "d2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

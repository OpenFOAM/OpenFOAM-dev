/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::fvPatch

Description
    A finiteVolume patch using a polyPatch and a fvBoundaryMesh

SourceFiles
    fvPatch.C
    fvPatchNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatch_H
#define fvPatch_H

#include "polyPatch.H"
#include "labelList.H"
#include "SubList.H"
#include "typeInfo.H"
#include "tmp.H"
#include "primitiveFields.H"
#include "SubField.H"
#include "fvPatchFieldsFwd.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvBoundaryMesh;
class objectRegistry;
class surfaceInterpolation;

/*---------------------------------------------------------------------------*\
                           Class fvPatch Declaration
\*---------------------------------------------------------------------------*/

class fvPatch
{
    // Private Data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;

        //- Reference to boundary mesh
        const fvBoundaryMesh& boundaryMesh_;


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

        //- Initialise the patches for moving points
        virtual void initMovePoints();

        //- Correct patches after moving points
        virtual void movePoints();


public:

    typedef fvBoundaryMesh BoundaryMesh;

    friend class fvBoundaryMesh;
    friend class surfaceInterpolation;

    //- Runtime type information
    TypeName(polyPatch::typeName_());


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvPatch,
            polyPatch,
            (const polyPatch& patch, const fvBoundaryMesh& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch and fvBoundaryMesh
        fvPatch(const polyPatch&, const fvBoundaryMesh&);

        //- Disallow default bitwise copy construction
        fvPatch(const fvPatch&);


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<fvPatch> New
        (
            const polyPatch&,
            const fvBoundaryMesh&
        );


    //- Destructor
    virtual ~fvPatch();


    // Member Functions

    // Access

            //- Return the polyPatch
            const polyPatch& patch() const
            {
                return polyPatch_;
            }

            //- Return name
            virtual const word& name() const
            {
                return polyPatch_.name();
            }

            //- Return start label of this patch in the polyMesh face list
            virtual label start() const
            {
                return polyPatch_.start();
            }

            //- Return size
            virtual label size() const
            {
                return polyPatch_.size();
            }

            //- Return true if this patch is coupled
            virtual bool coupled() const
            {
                return polyPatch_.coupled();
            }

            //- Return true if the given type is a constraint type
            static bool constraintType(const word& pt);

            //- Return a list of all the constraint patch types
            static wordList constraintTypes();

            //- Return the index of this patch in the fvBoundaryMesh
            label index() const
            {
                return polyPatch_.index();
            }

            //- Return boundaryMesh reference
            const fvBoundaryMesh& boundaryMesh() const
            {
                return boundaryMesh_;
            }

            //- Return the local object registry
            const objectRegistry& db() const;

            //- Slice list to patch
            template<class T>
            const typename List<T>::subList patchSlice(const List<T>& l) const
            {
                return typename List<T>::subList(l, size(), start());
            }

            //- Return faceCells
            virtual const labelUList& faceCells() const;


        // Access functions for geometrical data

            //- Return face centres
            const vectorField& Cf() const;

            //- Return neighbour cell centres
            tmp<vectorField> Cn() const;

            //- Return face area vectors
            const vectorField& Sf() const;

            //- Return face area magnitudes
            const scalarField& magSf() const;

            //- Return face normals
            tmp<vectorField> nf() const;

            //- Return cell-centre to face-centre vector
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre vector is returned
            virtual tmp<vectorField> delta() const;


        // Access functions for demand driven data

            //- Return patch weighting factors
            const scalarField& weights() const;

            //- Return the face - cell distance coeffient
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre distance coeffient is returned
            const scalarField& deltaCoeffs() const;


        // Evaluation functions

            //- Return given internal field next to patch as patch field
            template<class Type>
            tmp<Field<Type>> patchInternalField(const UList<Type>&) const;

            //- Return given internal field next to patch as patch field
            template<class Type>
            void patchInternalField(const UList<Type>&, Field<Type>&) const;

            //- Return the corresponding patchField of the named field
            template<class GeometricField, class Type>
            const typename GeometricField::Patch& patchField
            (
                const GeometricField&
            ) const;

            //- Return the corresponding patchField reference of the named field
            template<class GeometricField, class Type>
            typename GeometricField::Patch& patchField
            (
                GeometricField&
            ) const;

            //- Lookup and return the patchField of the named field from the
            //  local objectRegistry.
            template<class GeometricField, class Type>
            const typename GeometricField::Patch& lookupPatchField
            (
                const word& name
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

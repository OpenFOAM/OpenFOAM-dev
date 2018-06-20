/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::sampledSurface

Description
    An abstract class for surfaces with sampling.

    The constructors for the derived classes should generally start in a
    'expired' condition (ie, needsUpdate() == true) and rely on a
    subsequent call to the update() method to complete the initialization.
    Delaying the final construction as late as possible allows the
    construction of surfaces that may depend on intermediate calculation
    results (eg, iso-surfaces) and also avoids the unnecessary
    reconstruction of surfaces between sampling intervals.

    It is the responsibility of the caller to ensure that the surface
    update() is called before the surface is used.  The update() method
    implementation should do nothing when the surface is already
    up-to-date.

SourceFiles
    sampledSurface.C
    sampledSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurface_H
#define sampledSurface_H

#include "pointField.H"
#include "word.H"
#include "labelList.H"
#include "faceList.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceMesh.H"
#include "polyMesh.H"
#include "coordinateSystems.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class sampledSurface;

Ostream& operator<<(Ostream&, const sampledSurface&);


/*---------------------------------------------------------------------------*\
                      Class sampledSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledSurface
{
    // Private data

        //- Name of sample surface
        word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Do we intend to interpolate the information?
        const bool interpolate_;


    // Demand-driven data

        //- Face area vectors
        mutable vectorField* SfPtr_;

        //- Mag face area vectors
        mutable scalarField* magSfPtr_;

        //- Face centres
        mutable vectorField* CfPtr_;

        //- Total surface area
        mutable scalar area_;


    // Make geometric data

        //- Make Sf
        void makeSf() const;

        //- Make magSf
        void makeMagSf() const;

        //- Make Cf
        void makeCf() const;


    // Service methods

        //- Check field size matches surface size
        template<class Type>
        bool checkFieldSize(const Field<Type>&) const;

        //- Project field onto surface
        template<class ReturnType, class Type>
        void project
        (
            Field<ReturnType>&,
            const Field<Type>&
        ) const;

        //- Project field onto surface
        template<class ReturnType, class Type>
        void project
        (
            Field<ReturnType>&,
            const tmp<Field<Type>>&
        ) const;

        //- Project field onto surface
        template<class ReturnType, class Type>
        tmp<Field<ReturnType>> project(const tmp<Field<Type>>&) const;


protected:

    // Protected Member functions

        virtual void clearGeom() const;


public:

    //- Runtime type information
    TypeName("sampledSurface");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        sampledSurface,
        word,
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        ),
        (name, mesh, dict)
    );


    // iNew helper class

        //- Class used for the PtrLists read-construction
        class iNew
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<sampledSurface> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);

                return sampledSurface::New(name, mesh_, dict);
            }
        };


    // Constructors

        //- Construct from name, mesh
        sampledSurface
        (
            const word& name,
            const polyMesh&,
            const bool interpolate = false
        );

        //- Construct from dictionary
        sampledSurface
        (
            const word& name,
            const polyMesh&,
            const dictionary&
        );

        //- Clone
        autoPtr<sampledSurface> clone() const
        {
            NotImplemented;
            return autoPtr<sampledSurface>(nullptr);
        }


    // Selectors

        //- Return a reference to the selected surface
        static autoPtr<sampledSurface> New
        (
            const word& name,
            const polyMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~sampledSurface();


    // Member Functions

      // Access

        //- Access to the underlying mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Name of surface
        const word& name() const
        {
            return name_;
        }

        //- Interpolation requested for surface
        bool interpolate() const
        {
            return interpolate_;
        }

        //- Does the surface need an update?
        virtual bool needsUpdate() const = 0;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire() = 0;

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was required
        virtual bool update() = 0;

        //- Points of surface
        virtual const pointField& points() const = 0;

        //- Faces of surface
        virtual const faceList& faces() const = 0;

        //- Return face area vectors
        virtual const vectorField& Sf() const;

        //- Return face area magnitudes
        virtual const scalarField& magSf() const;

        //- Return face centres as vectorField
        virtual const vectorField& Cf() const;

        //- The total surface area
        scalar area() const;

        //- Integration of a field across the surface
        template<class Type>
        Type integrate(const Field<Type>&) const;

        //- Integration of a field across the surface
        template<class Type>
        Type integrate(const tmp<Field<Type>>&) const;

        //- Area-averaged value of a field across the surface
        template<class Type>
        Type average(const Field<Type>&) const;

        //- Area-averaged value of a field across the surface
        template<class Type>
        Type average(const tmp<Field<Type>>&) const;

        //- Project field onto surface
        tmp<Field<scalar>> project(const Field<scalar>&) const;

        //- Project field onto surface
        tmp<Field<scalar>> project(const Field<vector>&) const;

        //- Project field onto surface
        tmp<Field<vector>> project(const Field<sphericalTensor>&) const;

        //- Project field onto surface
        tmp<Field<vector>> project(const Field<symmTensor>&) const;

        //- Project field onto surface
        tmp<Field<vector>> project(const Field<tensor>&) const;

        //- Interpolate from points to cell centre
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> pointAverage
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pfld
        ) const;

        //- Sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const = 0;

        //- Surface sample field on surface
        virtual tmp<scalarField> sample
        (
            const surfaceScalarField&
        ) const;

        //- Surface Sample field on surface
        virtual tmp<vectorField> sample
        (
            const surfaceVectorField&
        ) const;

        //- Surface sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const surfaceSphericalTensorField&
        ) const;

        //- Surface sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const surfaceSymmTensorField&
        ) const;

        //- Surface sample field on surface
        virtual tmp<tensorField> sample
        (
            const surfaceTensorField&
        ) const;

        //- Interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const = 0;


        //- Interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const = 0;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const = 0;

        //- Interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const = 0;

        //- Interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const = 0;


      // Edit

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }


      // Write

        //- Write
        virtual void print(Ostream&) const;

        //- Ostream operator
        friend Ostream& operator<<(Ostream&, const sampledSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

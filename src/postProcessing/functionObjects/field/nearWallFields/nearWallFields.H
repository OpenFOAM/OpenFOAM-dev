/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::nearWallFields

Group
    grpFieldFunctionObjects

Description
    This function object samples near-patch volume fields

    Fields are stored
    - every time step the field is updated with new values
    - at output it writes the fields

    This functionObject can either be used
    - to calculate a new field as a  post-processing step or
    - since the fields are registered, used in another functionObject

    Example of function object specification:
    \verbatim
    nearWallFields1
    {
        type        nearWallFields;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        fields      ((p pNear)(U UNear));
        patches     (movingWall);
        distance    0.13;
    }
    \endverbatim

    \heading Function object usage
    \table
        Property | Description               | Required    | Default value
        type     | type name: nearWallFields | yes         |
        fields   | list of fields with correspoding output field names | yes |
        patches  | list of patches to sample | yes         |
        distance | distance from patch to sample | yes     |
    \endtable

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    nearWallFields.C
    IOnearWallFields.H

\*---------------------------------------------------------------------------*/

#ifndef nearWallFields_H
#define nearWallFields_H

#include "OFstream.H"
#include "volFields.H"
#include "Tuple2.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class nearWallFields Declaration
\*---------------------------------------------------------------------------*/

class nearWallFields
{
protected:

    // Protected data

        //- Name of this set of nearWallFields object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        // Read from dictionary

            //- Fields to process
            List<Tuple2<word, word> > fieldSet_;

            //- Patches to sample
            labelHashSet patchSet_;

            //- Distance away from wall
            scalar distance_;

            //- From original field to sampled result
            HashTable<word> fieldMap_;

            //- From resulting back to original field
            HashTable<word> reverseFieldMap_;


        // Calculated addressing

            //- From cell to seed patch faces
            labelListList cellToWalls_;

            //- From cell to tracked end point
            List<List<point> > cellToSamples_;

            //- Map from cell based data back to patch based data
            autoPtr<mapDistribute> getPatchDataMapPtr_;


        // Locally constructed fields

            PtrList<volScalarField> vsf_;
            PtrList<volVectorField> vvf_;
            PtrList<volSphericalTensorField> vSpheretf_;
            PtrList<volSymmTensorField> vSymmtf_;
            PtrList<volTensorField> vtf_;


    // Protected Member Functions

        //- Calculate addressing from cells back to patch faces
        void calcAddressing();

        template<class Type>
        void createFields
        (
            PtrList<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;

        //- Override boundary fields with sampled values
        template<class Type>
        void sampleBoundaryField
        (
            const interpolationCellPoint<Type>& interpolator,
            GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;

        template<class Type>
        void sampleFields
        (
            PtrList<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;


private:

        //- Disallow default bitwise copy construct
        nearWallFields(const nearWallFields&);

        //- Disallow default bitwise assignment
        void operator=(const nearWallFields&);

public:

    //- Runtime type information
    TypeName("nearWallFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        nearWallFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~nearWallFields();


    // Member Functions

        //- Return name of the nearWallFields object
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "nearWallFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

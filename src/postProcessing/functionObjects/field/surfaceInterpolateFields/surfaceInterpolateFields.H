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
    Foam::surfaceInterpolateFields

Group grpFieldFunctionObjects

Description
    This function object linearly interpolates volume fields to generate
    surface fields

    Fields are stored
    - every time step the field is updated with new values
    - at output it writes the fields

    This functionObject can either be used
    - to calculate a new field as a  post-processing step or
    - since the fields are registered, used in another functionObject

    Example of function object specification:
    \verbatim
    surfaceInterpolateFields1
    {
        type        surfaceInterpolateFields;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        fields      ((p pNear)(U UNear));
    }
    \endverbatim

    \heading Function object usage
    \table
        Property | Description               | Required    | Default value
        type     | type name: nearWallFields | yes         |
        fields   | list of fields with correspoding output field names | yes |
    \endtable


SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    surfaceInterpolateFields.C
    IOsurfaceInterpolateFields.H

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolateFields_H
#define surfaceInterpolateFields_H

#include "OFstream.H"
#include "surfaceFields.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                  Class surfaceInterpolateFields Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolateFields
{
protected:

    // Protected data

        //- Name of this set of surfaceInterpolateFields object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Fields to process
        //wordList fieldSet_;
        List<Tuple2<word, word> > fieldSet_;

        //- Locally constructed fields
        PtrList<surfaceScalarField> ssf_;
        PtrList<surfaceVectorField> svf_;
        PtrList<surfaceSphericalTensorField> sSpheretf_;
        PtrList<surfaceSymmTensorField> sSymmtf_;
        PtrList<surfaceTensorField> stf_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        surfaceInterpolateFields(const surfaceInterpolateFields&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceInterpolateFields&);

        template<class Type>
        void interpolateFields
        (
            PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> >&
        ) const;


public:

    //- Runtime type information
    TypeName("surfaceInterpolateFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        surfaceInterpolateFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~surfaceInterpolateFields();


    // Member Functions

        //- Return name of the surfaceInterpolateFields object
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
#   include "surfaceInterpolateFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

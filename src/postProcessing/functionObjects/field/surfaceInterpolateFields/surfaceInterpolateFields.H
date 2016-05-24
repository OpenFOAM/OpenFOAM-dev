/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
        libs ("libfieldFunctionObjects.so");
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
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    surfaceInterpolateFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceInterpolateFields_H
#define functionObjects_surfaceInterpolateFields_H

#include "fvMeshFunctionObject.H"
#include "surfaceFieldsFwd.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class surfaceInterpolateFields Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolateFields
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Fields to process
        List<Tuple2<word, word>> fieldSet_;

        //- Locally constructed fields
        PtrList<surfaceScalarField> ssf_;
        PtrList<surfaceVectorField> svf_;
        PtrList<surfaceSphericalTensorField> sSpheretf_;
        PtrList<surfaceSymmTensorField> sSymmtf_;
        PtrList<surfaceTensorField> stf_;


    // Protected Member Functions

        template<class Type>
        void interpolateFields
        (
            PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>&
        ) const;


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        surfaceInterpolateFields(const surfaceInterpolateFields&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceInterpolateFields&);


public:

    //- Runtime type information
    TypeName("surfaceInterpolateFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        surfaceInterpolateFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceInterpolateFields();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the interpolated fields
        virtual bool execute(const bool postProcess = false);

        //- Write the interpolated fields
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceInterpolateFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::turbulenceFields

Group
    grpFieldFunctionObjects

Description
    This function object stores turbulence fields on the mesh database for
    further manipulation.

    Fields are stored as copies of the original, with the prefix
    "tubulenceModel:", e.g.

        turbulenceModel:R

    Example of function object specification:
    \verbatim
    turbulenceFields1
    {
        type        turbulenceFields;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
        fields
        (
            R
            devRhoReff
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: processorField | yes       |
        fields       | fields to store (see below) | yes     |
    \endtable

    Where \c fields can include:
    \plaintable
        R           | Stress tensor
        devRhoReff  |
        mut         | turbulence viscosity (compressible)
        muEff       | effective turbulence viscosity (compressible)
        alphat      | turbulence thermal diffusivity (compressible)
        alphaEff    | effective turbulence thermal diffusivity (compressible)
        devReff     |
        nut         | turbulence viscosity (incompressible)
        nuEff       | effective turbulence viscosity (incompressible)
    \endplaintable

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    turbulenceFields.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceFields_H
#define turbulenceFields_H

#include "HashSet.H"
#include "IOobject.H"
#include "NamedEnum.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                      Class turbulenceFields Declaration
\*---------------------------------------------------------------------------*/

class turbulenceFields
{
public:

    enum compressibleField
    {
        cfR,
        cfDevRhoReff,
        cfMut,
        cfMuEff,
        cfAlphat,
        cfAlphaEff
    };
    static const NamedEnum<compressibleField, 6> compressibleFieldNames_;

    enum incompressibleField
    {
        ifR,
        ifDevReff,
        ifNut,
        ifNuEff
    };
    static const NamedEnum<incompressibleField, 4> incompressibleFieldNames_;

    static const word modelName;


protected:

    // Protected data

        //- Name of this set of turbulenceFields object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Fields to load
        wordHashSet fieldSet_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        turbulenceFields(const turbulenceFields&);

        //- Disallow default bitwise assignment
        void operator=(const turbulenceFields&);

        //- Return true if compressible turbulence model is identified
        bool compressible();

        //- Process the turbulence field
        template<class Type>
        void processField
        (
            const word& fieldName,
            const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvalue
        );


public:

    //- Runtime type information
    TypeName("turbulenceFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        turbulenceFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~turbulenceFields();


    // Member Functions

        //- Return name of the turbulenceFields object
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
    #include "turbulenceFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

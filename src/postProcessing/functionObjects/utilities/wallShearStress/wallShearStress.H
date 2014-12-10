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
    Foam::wallShearStress

Group
    grpUtilitiesFunctionObjects

Description
    This function object evaluates and outputs the shear stress at wall
    patches.  The result is written as a volVectorField to time directories as
    field 'wallShearStress'

        \f[
            Stress = R \dot n
        \f]

    where
    \vartable
        R       | stress tensor
        n       | patch normal vector (into the domain)
    \endvartable

    The shear stress (symmetrical) tensor field is retrieved from the
    turbulence model.  All wall patches are included by default; to restrict
    the calculation to certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    wallShearStress1
    {
        type        wallShearStress;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
        patches     (".*Wall");
    }
    \endverbatim

    \heading Function object usage
    \table
        Property | Description               | Required    | Default value
        type     | type name: wallShearStress | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

SourceFiles
    wallShearStress.C
    IOwallShearStress.H

\*---------------------------------------------------------------------------*/

#ifndef wallShearStress_H
#define wallShearStress_H

#include "functionObjectFile.H"
#include "volFieldsFwd.H"
#include "Switch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class wallShearStress Declaration
\*---------------------------------------------------------------------------*/

class wallShearStress
:
    public functionObjectFile
{
protected:

    // Protected data

        //- Name of this set of wallShearStress object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the shear stress
        void calcShearStress
        (
            const fvMesh& mesh,
            const volSymmTensorField& Reff,
            volVectorField& shearStress
        );

        //- Disallow default bitwise copy construct
        wallShearStress(const wallShearStress&);

        //- Disallow default bitwise assignment
        void operator=(const wallShearStress&);


public:

    //- Runtime type information
    TypeName("wallShearStress");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        wallShearStress
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~wallShearStress();


    // Member Functions

        //- Return name of the set of wallShearStress
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the wallShearStress data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the wallShearStress and write
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

#endif

// ************************************************************************* //

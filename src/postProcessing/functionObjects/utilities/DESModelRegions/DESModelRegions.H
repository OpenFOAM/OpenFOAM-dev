/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::DESModelRegions

Group
    grpUtilitiesFunctionObjects

Description
    This function object writes out an indicator field for DES turbulence
    calculations, that is:

    - 0 for RAS regions
    - 1 for LES regions

SourceFiles
    DESModelRegions.C
    IODESModelRegions.H

\*---------------------------------------------------------------------------*/

#ifndef DESModelRegions_H
#define DESModelRegions_H

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
                       Class DESModelRegions Declaration
\*---------------------------------------------------------------------------*/

class DESModelRegions
:
    public functionObjectFile
{
protected:

    // Protected data

        //- Name of this set of DESModelRegions object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Disallow default bitwise copy construct
        DESModelRegions(const DESModelRegions&);

        //- Disallow default bitwise assignment
        void operator=(const DESModelRegions&);


public:

    //- Runtime type information
    TypeName("DESModelRegions");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        DESModelRegions
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~DESModelRegions();


    // Member Functions

        //- Return name of the set of DESModelRegions
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the DESModelRegions data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the DESModelRegions and write
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

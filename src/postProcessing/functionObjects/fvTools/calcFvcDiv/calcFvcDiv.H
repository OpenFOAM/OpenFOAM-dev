/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::calcFvcDiv

Group
    grpFVFunctionObjects

Description
    This function object calculates the divergence of a field.  The operation is
    limited to surfaceScalarFields and volumeVector fields, and the output is a
    volume scalar field.

SourceFiles
    calcFvcDiv.C
    IOcalcFvcDiv.H

\*---------------------------------------------------------------------------*/

#ifndef calcFvcDiv_H
#define calcFvcDiv_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldFwd.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class dimensionSet;

/*---------------------------------------------------------------------------*\
                         Class calcFvcDiv Declaration
\*---------------------------------------------------------------------------*/

class calcFvcDiv
{
    // Private data

        //- Name of this calcFvcDiv object
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of field to process
        word fieldName_;

        //- Name of result field
        word resultName_;


    // Private Member Functions

        //- Helper function to create/store/return the divergence field
        volScalarField& divField
        (
            const word& gradName,
            const dimensionSet& dims
        );

        //- Helper function to calculate the divergence of different field types
        template<class FieldType>
        void calcDiv
        (
            const word& fieldName,
            const word& resultName,
            bool& processed
        );

        //- Disallow default bitwise copy construct
        calcFvcDiv(const calcFvcDiv&);

        //- Disallow default bitwise assignment
        void operator=(const calcFvcDiv&);


public:

    //- Runtime type information
    TypeName("calcFvcDiv");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        calcFvcDiv
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~calcFvcDiv();


    // Member Functions

        //- Return name of the set of calcFvcDiv
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the calcFvcDiv data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the calcFvcDiv and write
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
    #include "calcFvcDivTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

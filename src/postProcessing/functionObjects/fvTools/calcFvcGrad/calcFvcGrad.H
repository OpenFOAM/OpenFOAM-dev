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
    Foam::calcFvcGrad

Group
    grpFVFunctionObjects

Description
    This function object calculates the gradient of a field.  The operation is
    limited to scalar and vector volume or surface fields, and the output is a
    volume vector or tensor field.

SourceFiles
    calcFvcGrad.C
    IOcalcFvcGrad.H

\*---------------------------------------------------------------------------*/

#ifndef calcFvcGrad_H
#define calcFvcGrad_H

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
                         Class calcFvcGrad Declaration
\*---------------------------------------------------------------------------*/

class calcFvcGrad
{
    // Private data

        //- Name of this calcFvcGrad object
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

        //- Helper function to create/store/return the gradient field
        template<class Type>
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >&
        gradField(const word& gradName, const dimensionSet& dims);

        //- Helper function to calculate the gradient of different field types
        template<class Type>
        void calcGrad
        (
            const word& fieldName,
            const word& resultName,
            bool& processed
        );

        //- Disallow default bitwise copy construct
        calcFvcGrad(const calcFvcGrad&);

        //- Disallow default bitwise assignment
        void operator=(const calcFvcGrad&);


public:

    //- Runtime type information
    TypeName("calcFvcGrad");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        calcFvcGrad
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~calcFvcGrad();


    // Member Functions

        //- Return name of the set of calcFvcGrad
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the calcFvcGrad data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the calcFvcGrad and write
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
    #include "calcFvcGradTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

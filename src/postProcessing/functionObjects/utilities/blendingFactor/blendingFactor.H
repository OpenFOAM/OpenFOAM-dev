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
    Foam::blendingFactor

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs the blendingFactor as used by
    the bended convection schemes.  The output is a volume field (cells) whose
    value is calculated via the maximum blending factor for any cell face.


SourceFiles
    blendingFactor.C
    IOblendingFactor.H

\*---------------------------------------------------------------------------*/

#ifndef blendingFactor_H
#define blendingFactor_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
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

/*---------------------------------------------------------------------------*\
                          Class blendingFactor Declaration
\*---------------------------------------------------------------------------*/

class blendingFactor
{
    // Private data

        //- Name of this set of blendingFactor objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Field name
        word fieldName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        blendingFactor(const blendingFactor&);

        //- Disallow default bitwise assignment
        void operator=(const blendingFactor&);

        //- Return the blending factor field from the database
        template<class Type>
        volScalarField& factor
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        );

        //- Calculate the blending factor
        template<class Type>
        void calc();


public:

    //- Runtime type information
    TypeName("blendingFactor");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        blendingFactor
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~blendingFactor();


    // Member Functions

        //- Return name of the set of blendingFactor
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the blendingFactor data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the blendingFactor and write
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
    #include "blendingFactorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

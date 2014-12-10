/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::regionModelFunctionObject

Description
    Region model function object base class

SourceFiles
    regionModelFunctionObject.C
    regionModelFunctionObjectNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionModelFunctionObject_H
#define regionModelFunctionObject_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regionModels
{

class regionModel;

/*---------------------------------------------------------------------------*\
                  Class regionModelFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class regionModelFunctionObject
{
protected:

    // Protected data

        //- Dictionary
        dictionary dict_;

        //- Reference to the region model
        regionModel& owner_;

        //- Model type name
        word modelType_;


public:

    //- Runtime type information
    TypeName("regionModelFunctionObject");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        regionModelFunctionObject,
        dictionary,
        (
            const dictionary& dict,
            regionModel& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        regionModelFunctionObject(regionModel& owner);

        //- Construct from dictionary
        regionModelFunctionObject
        (
            const dictionary& dict,
            regionModel& owner,
            const word& modelType
        );

        //- Construct copy
        regionModelFunctionObject(const regionModelFunctionObject& ppm);

        //- Construct and return a clone
        virtual autoPtr<regionModelFunctionObject> clone() const
        {
            return autoPtr<regionModelFunctionObject>
            (
                new regionModelFunctionObject(*this)
            );
        }


    //- Destructor
    virtual ~regionModelFunctionObject();


    //- Selector
    static autoPtr<regionModelFunctionObject> New
    (
        const dictionary& dict,
        regionModel& owner,
        const word& modelType
    );


    // Member Functions

        // Evaluation

            //- Pre-evolve region hook
            virtual void preEvolveRegion();

            //- Post-evolve region hook
            virtual void postEvolveRegion();

    // I-O

        //- write
        virtual void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::CourantNo

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs the Courant number as a
    volScalarField.  The field is stored on the mesh database so that it can
    be retrieved and used for other applications.

SourceFiles
    CourantNo.C
    IOCourantNo.H

\*---------------------------------------------------------------------------*/

#ifndef CourantNo_H
#define CourantNo_H

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
                          Class CourantNo Declaration
\*---------------------------------------------------------------------------*/

class CourantNo
{
    // Private data

        //- Name of this set of CourantNo objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;


    // Private Member Functions

        //- Return 1 if the flux field is volumetric, otherwise return rho
        //  from the database
        tmp<volScalarField> rho(const surfaceScalarField& p) const;

        //- Disallow default bitwise copy construct
        CourantNo(const CourantNo&);

        //- Disallow default bitwise assignment
        void operator=(const CourantNo&);


public:

    //- Runtime type information
    TypeName("CourantNo");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        CourantNo
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~CourantNo();


    // Member Functions

        //- Return name of the set of CourantNo
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the CourantNo data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the CourantNo and write
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

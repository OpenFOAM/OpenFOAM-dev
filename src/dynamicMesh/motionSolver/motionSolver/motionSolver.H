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
    Foam::motionSolver

Description
    Virtual base class for mesh motion solver.

SourceFiles
    motionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef motionSolver_H
#define motionSolver_H

#include "IOdictionary.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class motionSolver Declaration
\*---------------------------------------------------------------------------*/

class motionSolver
:
    public IOdictionary
{

private:

    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Model coefficients dictionary
        dictionary coeffDict_;


    // Private Member Functions

        //- De-register object if registered and assign to current
        static IOobject stealRegistration(const IOdictionary& dict);

public:

    //- Runtime type information
    TypeName("motionSolver");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            motionSolver,
            dictionary,
            (const polyMesh& mesh, const IOdictionary& dict),
            (mesh, dict)
        );


    // Selectors

        //- Select constructed from polyMesh
        static autoPtr<motionSolver> New(const polyMesh&);

        //- Select constructed from polyMesh and dictionary. If dictionary
        //  was registered this will 'steal' that registration.
        static autoPtr<motionSolver> New(const polyMesh&, const IOdictionary&);


    // Constructors

        //- Construct from polyMesh
        motionSolver(const polyMesh& mesh);

        //- Construct from polyMesh and dictionary and type.
        motionSolver
        (
            const polyMesh& mesh,
            const IOdictionary&,
            const word& type
        );


    //- Destructor
    virtual ~motionSolver();


    // Member Functions

        //- Return reference to mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Const access to the coefficients dictionary
        const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Provide new points for motion.  Solves for motion
        virtual tmp<pointField> newPoints();

        //- Provide current points for motion.  Uses current motion field
        virtual tmp<pointField> curPoints() const = 0;

        virtual void twoDCorrectPoints(pointField&) const;

        //- Solve for motion
        virtual void solve() = 0;

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&) = 0;

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&) = 0;

        //- Write state using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;

        //- Read dynamicMeshDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

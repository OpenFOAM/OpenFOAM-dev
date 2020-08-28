/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "pointField.H"
#include "Time.H"
#include "polyMesh.H"
#include "dictionaryEntry.H"

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
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Model coefficients dictionary
        dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("motionSolver");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            motionSolver,
            dictionary,
            (const polyMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Selectors

        //- Select constructed from polyMesh and dictionary
        static autoPtr<motionSolver> New(const polyMesh&, const dictionary&);

        //- Class used for the construction of PtrLists of motionSolvers
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh);

            autoPtr<motionSolver> operator()(Istream& is) const;
        };


    // Constructors

        //- Construct from polyMesh and dictionary and type.
        motionSolver
        (
            const polyMesh& mesh,
            const dictionary&,
            const word& type
        );

        //- Clone function
        virtual autoPtr<motionSolver> clone() const;


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

        //- Optionally write motion state information for restart
        virtual bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

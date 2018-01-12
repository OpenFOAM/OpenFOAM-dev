/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::convergenceControl

Description
    Convergence control class. Provides methods to check the convergence of the
    time loop against an absolute residual tolerance.

SourceFiles
    convergenceControl.C
    convergenceControlTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef convergenceControl_H
#define convergenceControl_H

#include "fvMesh.H"
#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class convergenceControl Declaration
\*---------------------------------------------------------------------------*/

class convergenceControl
{
public:

    //- Classes

        //- Residual data structure
        struct residualData
        {
            wordRe name;
            scalar absTol;
        };


    // Static Functions

        //- Get the initial residuals for the first and the i-th solves in this
        //  time-step
        static void getInitialResiduals
        (
            const fvMesh& mesh,
            const word& fieldName,
            const label solvei,
            ITstream& data,
            scalar& r0,
            scalar& r
        );


    // Static Template Functions

        //- Return the index of the named field in residual control data, or -1
        //  if not present
        template<class ResidualData>
        static label residualControlIndex
        (
            const word& fieldName,
            const List<ResidualData>& residualControl,
            const bool useRegEx=true
        );

        //- Get the initial residuals for the first and the i-th solves in this
        //  time-step
        template<class Type>
        static void getInitialTypeResiduals
        (
            const fvMesh& mesh,
            const word& fieldName,
            const label solvei,
            ITstream& data,
            scalar& r0,
            scalar& r
        );



private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        convergenceControl(const convergenceControl&);

        //- Disallow default bitwise assignment
        void operator=(const convergenceControl&);


protected:

    // Protected data

        //- Reference to the solution control
        const solutionControl& control_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("convergenceControl");


    // Constructors

        //- Construct from a solution control
        convergenceControl(const solutionControl& control);


    //- Destructor
    virtual ~convergenceControl();


    // Member Functions

        // Evolution

            //- Return true if residual controls are present
            virtual bool hasResidualControls() const = 0;

            //- Return true if all convergence checks are satisfied
            virtual bool criteriaSatisfied() const = 0;

            //- Flag to indicate whether convergance has been reached
            bool converged();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "convergenceControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

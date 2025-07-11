/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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
    Foam::functionObjects::age

Description
    Calculates and writes out the time taken for a particle to travel from an
    inlet to the location. Solves the following equation when incompressible:
    \f[
        \div (\phi t) = 1
    \f]
    where:
    \vartable
        t    | Age [s]
        \phi | Volumetric flux [m^3/s]
    \endvartable
    Boundary conditions are generated automatically as zeroGradient on all
    walls and inletOutlet everywhere else.

Usage
    \table
        Property  | Description                   | Required | Default value
        phi       | The name of the flux field    | no       | phi
        rho       | The name of the density field | no       | rho
        nCorr     | The maximum number of correctors | no       | 5
        schemesField | The name of the field from which schemes are taken | \\
        no | age
        diffusion | Switch to turn on/off the diffusion term | no | off
        tolerance  | Solver residual control       | no       | 1e-5
    \endtable

    \verbatim
    age1
    {
        type            age;
        libs            ("libsolverFunctionObjects.so");

        executeControl  writeTime;
        writeControl    writeTime;

        schemesField    k;
    }
    \endverbatim

SourceFiles
    age.C

\*---------------------------------------------------------------------------*/

#ifndef age_functionObject_H
#define age_functionObject_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                             Class age Declaration
\*---------------------------------------------------------------------------*/

class age
:
    public fvMeshFunctionObject
{
    // Private Data

        //- The name of the flux field
        word phiName_;

        //- The name of the density field
        word rhoName_;

        //- Number of corrections
        int nCorr_;

        //- Name of field from which schemes are taken
        word schemesField_;

        //- Switch to turn on/off the diffusion term
        Switch diffusion_;

        //- Convergence tolerance
        scalar tolerance_;


    // Private Member Functions

        //- The list of patch types for the age field
        wordList patchTypes() const;

        //- Check convergence
        bool converged(const int nCorr, const scalar initialResidual) const;


public:

    //- Runtime type information
    TypeName("age");


    // Constructors

        //- Construct from Time and dictionary
        age
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~age();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

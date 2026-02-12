/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::functionObjects::sectionalForceProbes

Description
    This function calculates the fluid forces and moments at a number of
    section-planes through a given set of patches. It writes the result into a
    log file.

    The cut-planes are defined by a single normal and a number of points.
    Forces and moments are calculated on the side of the planes to which the
    normal points. If the "other" side also needs to be generated then the
    function can be executed with the normal reversed.

    Note that this function only generates the contribution from the fluid to
    the sectional forces and moments in the object. To compute the complete
    sectional force would require supporting and/or dynamic forces to be
    included also.

    Example of function object specification:
    \verbatim
    sectionalForceProbes1
    {
        type            sectionalForceProbes;
        libs            ("libforces.so");
        patches         (hull);
        points
        (
            (-0.8 0 0)
            ( 0   0 0)
            ( 0.8 0 0)
        );
        normal          (1 0 0);
    }
    \endverbatim

Usage
    \table
        Property     | Description              | Required   | Default value
        type         | Type name: sectionalForceProbes | yes |
        patches      | Patches included in the forces calculation | yes |
        points       | Points in the cut-planes | yes        |
        normal       | Normal to the cut-planes | yes        |
        p            | Pressure field name      | no         | p
        U            | Velocity field name      | no         | U
        rho          | Density field name       | no         | rho
        phase        | Phase name for phase-fraction | no    |
        rhoRef       | Reference density        | if incompressible |
        pRef         | Reference pressure       | no         | 0
    \endtable

SourceFiles
    sectionalForceProbes.C

\*---------------------------------------------------------------------------*/

#ifndef sectionalForceProbes_functionObject_H
#define sectionalForceProbes_functionObject_H

#include "sectionalForcesBase.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class sectionalForceProbes Declaration
\*---------------------------------------------------------------------------*/

class sectionalForceProbes
:
    public sectionalForcesBase,
    public logFiles
{
    // Private Data

        //- The normal to the selection planes
        vector normal_;

        //- Points in the section planes
        List<point> points_;

        //- The order the points when sorted by distance in the normal
        labelList pointOrders_;


protected:

    // Protected Member Functions

        //- Return the normal to the section planes
        virtual vector normal() const;

        //- Return the origin of the plot
        virtual point origin() const;

        //- Return the distances from the origin to the section planes
        virtual tmp<scalarField> distances() const;

        //- Output file header information
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("sectionalForceProbes");


    // Constructors

        //- Construct from Time and dictionary
        sectionalForceProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sectionalForceProbes(const sectionalForceProbes&) = delete;


    //- Destructor
    virtual ~sectionalForceProbes();


    // Member Functions

        //- Read the sectionalForceProbes data
        virtual bool read(const dictionary&);

        //- Write the sectionalForceProbes
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sectionalForceProbes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

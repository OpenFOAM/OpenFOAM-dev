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
    Foam::functionObjects::sectionalForceGraph

Description
    This function calculates the fluid forces and moments at a series of
    section-planes through a given set of patches. It writes the result into a
    graph which can be used to construct part of a shear force and bending
    moment diagram. Points in the graph are graded to approximately match the
    resolution of the mesh.

    The cut-planes are defined by a normal, and the line of the graph is
    defined by an origin and extends in the direction of the normal. Forces and
    moments are calculated on the side of the planes to which the normal
    points. If the "other" side also needs to be generated then the function
    can be executed with the normal reversed.

    Note that to construct full shear force and bending moment diagrams would
    require supporting and/or dynamic forces to be included also. This function
    only generates the contribution from the surrounding fluid.

    Example of function object specification:
    \verbatim
    sectionalForceGraph1
    {
        type            sectionalForceGraph;
        libs            ("libforces.so");
        patches         (hull);
        origin          (0 0 0);
        normal          (1 0 0);
        nPoints         20;
        axis            x;
        setFormat       raw;
    }
    \endverbatim

Usage
    \table
        Property     | Description              | Required   | Default value
        type         | Type name: sectionalForceGraph | yes  |
        patches      | Patches included in the forces calculation | yes |
        origin       | Origin of the plot       | yes        |
        normal       | Normal to the section planes | yes    |
        nPoints      | Number of points in the plot | yes    |
        axis         | Component of the position to plot against | yes |
        setFormat    | Graph output format      | yes        |
        p            | Pressure field name      | no         | p
        U            | Velocity field name      | no         | U
        rho          | Density field name       | no         | rho
        phase        | Phase name for phase-fraction | no    |
        rhoRef       | Reference density        | if incompressible |
        pRef         | Reference pressure       | no         | 0
    \endtable

SourceFiles
    sectionalForceGraph.C

\*---------------------------------------------------------------------------*/

#ifndef sectionalForceGraph_functionObject_H
#define sectionalForceGraph_functionObject_H

#include "sectionalForcesBase.H"
#include "coordSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class sectionalForceGraph Declaration
\*---------------------------------------------------------------------------*/

class sectionalForceGraph
:
    public sectionalForcesBase
{
    // Private Data

        //- The normal to the selection planes
        vector normal_;

        //- The origin of the plot
        vector origin_;

        //- The distances from the origin to the section planes
        mutable autoPtr<scalarField> distancesPtr_;

        //- The number of points in the plot
        label nPoints_;

        //- The direction over which to plot the results
        coordSet::axisType axis_;

        //- File writer
        autoPtr<setWriter> formatter_;

        //- Number of optimisation iterations. Default is 2.
        label nOptimiseIter_;


    // Private Member Functions

        //- Clear the cached weights
        virtual void clear();


protected:

    // Protected Member Functions

        //- Return the normal to the section planes
        virtual vector normal() const;

        //- Return the origin of the plot
        virtual point origin() const;

        //- Return the distances from the origin to the section planes
        virtual tmp<scalarField> distances() const;


public:

    //- Runtime type information
    TypeName("sectionalForceGraph");


    // Constructors

        //- Construct from Time and dictionary
        sectionalForceGraph
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sectionalForceGraph(const sectionalForceGraph&) = delete;


    //- Destructor
    virtual ~sectionalForceGraph();


    // Member Functions

        //- Read the sectionalForceGraph data
        virtual bool read(const dictionary&);

        //- Write the sectionalForceGraph
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sectionalForceGraph&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

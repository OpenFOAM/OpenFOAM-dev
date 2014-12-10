/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::fv::explicitPorositySource

Description
    Explicit porosity source

    \heading Source usage

    Example usage, here employing the Darcy-Forchheimer model:
    \verbatim
    explicitPorositySourceCoeffs
    {
        type            DarcyForchheimer;

        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (5e7 -1000 -1000);
            f   f [0 -1 0 0 0 0 0] (0 0 0);

            coordinateSystem
            {
                e1  (0.70710678 0.70710678 0);
                e2  (0 0 1);
            }
        }
    }
    \endverbatim

Note:
    The porous region must be selected as a cellZone.

SourceFiles
    explicitPorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef explicitPorositySource_H
#define explicitPorositySource_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class explicitPorositySource Declaration
\*---------------------------------------------------------------------------*/

class explicitPorositySource
:
    public option
{

protected:

    // Protected data

        //- Run-time selectable porosity model
        autoPtr<porosityModel> porosityPtr_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        explicitPorositySource(const explicitPorositySource&);

        //- Disallow default bitwise assignment
        void operator=(const explicitPorositySource&);


public:

    //- Runtime type information
    TypeName("explicitPorositySource");


    // Constructors

        //- Construct from components
        explicitPorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~explicitPorositySource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

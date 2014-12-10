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
    Foam::porosityModelList

Description
    List container for porosity models

SourceFiles
    porosityModelList.C

\*---------------------------------------------------------------------------*/

#ifndef porosityModelList_H
#define porosityModelList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "porosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class porosityModelList;
Ostream& operator<<(Ostream& os, const porosityModelList& models);

/*---------------------------------------------------------------------------*\
                      Class porosityModelList Declaration
\*---------------------------------------------------------------------------*/

class porosityModelList
:
    public PtrList<porosityModel>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        porosityModelList(const porosityModelList&);

        //- Disallow default bitwise assignment
        void operator=(const porosityModelList&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;


public:

    //- Constructor
    porosityModelList(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    ~porosityModelList();


    // Member Functions

        //- Return active status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Add resistance
        void addResistance(fvVectorMatrix& UEqn);

        //- Add resistance
        void addResistance
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC = true
        );


        // I-O

            //- Read dictionary
            bool read(const dictionary& dict);

            //- Write data to Ostream
            bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const porosityModelList& models
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

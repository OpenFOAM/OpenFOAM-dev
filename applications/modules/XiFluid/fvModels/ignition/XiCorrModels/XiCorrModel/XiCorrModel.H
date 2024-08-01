/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::XiCorrModel

Description
    Base class for ignition kernel flame wrinkling \c Xi correction

Usage
    Example usage:
    \verbatim
    XiCorr
    {
        type            <type>;
        select          all;
        bMin            0.05;
    }
    \endverbatim

    Where:
    \table
        Property         | Description              | Required  | Default value
        select           | Correction cell set selection | yes  |
        bMin             | Min b below which no correction | no | 0.01
    \endtable

SourceFiles
    XiCorrModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiCorrModel_H
#define XiCorrModel_H

#include "fvCellSet.H"
#include "volFieldsFwd.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiCorrModel Declaration
\*---------------------------------------------------------------------------*/

class XiCorrModel
{
    // Private Data

        //- The set of Xi correction cells
        fvCellSet set_;

        //- Minimum b below which the Xi correction is not applied
        //  Defaults to 0.01
        dimensionedScalar bMin_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("XiCorrModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiCorrModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from mesh and dictionary
        XiCorrModel(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        XiCorrModel(const XiCorrModel&) = delete;


    // Selectors

        //- Return a reference to the selected XiCorrModel model
        static autoPtr<XiCorrModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~XiCorrModel();


    // Member Functions

        //- Return the area of the ignition kernel
        //  calculated from the given kernel volume
        virtual dimensionedScalar Ak(const dimensionedScalar& Vk) const = 0;

        virtual void XiCorr
        (
            volScalarField& Xi,
            const volScalarField& b,
            const volScalarField& mgb
        ) const;


        // Mesh motion

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // IO

        //- Update properties from the given XiProperties dictionary
        bool read(const dictionary& XiProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const XiCorrModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

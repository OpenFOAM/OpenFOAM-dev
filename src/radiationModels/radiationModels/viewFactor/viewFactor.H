/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::radiationModels::viewFactor

Description
    View factor radiation model. The system solved is: C q = b
    where:
            Cij  = deltaij/Ej - (1/Ej - 1)Fij
            q    = heat flux
            b    = A eb - Ho
    and:
            eb   = sigma*T^4
            Ej   = emissivity
            Aij  = deltaij - Fij
            Fij  = view factor matrix


SourceFiles
    viewFactor.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_viewFactor_H
#define radiation_viewFactor_H

#include "radiationModel.H"
#include "singleCellFvMesh.H"
#include "scalarMatrices.H"
#include "globalIndex.H"
#include "scalarListIOList.H"
#include "distributionMap.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                         Class viewFactor Declaration
\*---------------------------------------------------------------------------*/

class viewFactor
:
    public radiationModel
{
    // Private Data

        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Map distributed
        autoPtr<distributionMap> map_;

        //- Coarse mesh
        singleCellFvMesh coarseMesh_;

        //- Net radiative heat flux [W/m^2]
        volScalarField qr_;

        //- View factor matrix
        autoPtr<scalarSquareMatrix> Fmatrix_;

        //- Inverse of C matrix
        autoPtr<scalarSquareMatrix> CLU_;

        //- Selected patches
        labelList selectedPatches_;

        //- Total global coarse faces
        label totalNCoarseFaces_;

        //- Total local coarse faces
        label nLocalCoarseFaces_;

        //- Constant emissivity
        bool constEmissivity_;

        //- Iterations Counter
        label iterCounter_;

        //- Pivot Indices for LU decomposition
        labelList pivotIndices_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Insert view factors into main matrix
        void insertMatrixElements
        (
            const globalIndex& index,
            const label fromProci,
            const labelListList& globalFaceFaces,
            const scalarListList& viewFactors,
            scalarSquareMatrix& matrix
        );


public:

    //- Runtime type information
    TypeName("viewFactor");


    // Constructors

        //- Construct from components
        viewFactor(const volScalarField& T);

        //- Construct from components
        viewFactor(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        viewFactor(const viewFactor&) = delete;


    //- Destructor
    virtual ~viewFactor();


    // Member Functions

        // Edit

            //- Solve system of equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;


    // Access

        //- Const access to total radiative heat flux field
        inline const volScalarField& qr() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viewFactor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "viewFactorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

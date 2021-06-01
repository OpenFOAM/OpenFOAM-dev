/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::chemistryTabulationMethods::ISAT

Description
    Implementation of the ISAT (In-situ adaptive tabulation), for chemistry
    calculation.

    Reference:
    \verbatim
        Pope, S. B. (1997).
        Computationally efficient implementation of combustion chemistry using
        in situ adaptive tabulation.
        Combustion Theory and Modelling, 1, 41-63.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ISAT_H
#define ISAT_H

#include "binaryTree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryTabulationMethods
{

/*---------------------------------------------------------------------------*\
                            Class ISAT Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class ISAT
:
    public chemistryTabulationMethod<ThermoType>
{
    // Private Data

        //- List of the stored 'points' organised in a binary tree
        binaryTree<ThermoType> chemisTree_;

        //- List of scale factors for species, temperature and pressure
        scalarField scaleFactor_;

        const Time& runTime_;

        //- Lifetime (number of time steps) of a stored point
        label chPMaxLifeTime_;

        //- Maximum number of growths before removing from the tree
        label maxGrowth_;

        //- Check the binary tree for leafs to remove every interval
        label checkEntireTreeInterval_;

        //- Factor that multiply the ideal depth of a binary tree to decide
        // whether to try to balance of not
        scalar maxDepthFactor_;

        //- Minimal size before trying to balance the tree
        label minBalanceThreshold_;

        //- After a failed primary retrieve, look in the MRU list
        Switch MRURetrieve_;

        //- Most Recently Used (MRU) list of chemPoint
        SLList<chemPointISAT<ThermoType>*> MRUList_;

        //- Maximum size of the MRU list
        label maxMRUSize_;

        //- Store a pointer to the last chemPointISAT found
        chemPointISAT<ThermoType>* lastSearch_;

        //- Switch to allow growth (on by default)
        Switch growPoints_;

        // Statistics on ISAT usage
        label nRetrieved_;
        label nGrowth_;
        label nAdd_;

        autoPtr<OFstream> nRetrievedFile_;
        autoPtr<OFstream> nGrowthFile_;
        autoPtr<OFstream> nAddFile_;
        autoPtr<OFstream> sizeFile_;

        bool cleaningRequired_;


    // Private Member Functions

        //- Add a chemPoint to the MRU list
        void addToMRU(chemPointISAT<ThermoType>* phi0);

        //- Compute and return the mapping of the composition phiq
        //  Input : phi0 the nearest chemPoint used in the linear interpolation
        //  phiq the composition of the query point for which we want to
        //  compute the mapping
        //  Rphiq the mapping of the new composition point (given as empty)
        //  Output: void (the mapping is stored in the Rphiq array)
        //  Rphiq = Rphi0 + A * (phiq-phi0)
        void calcNewC
        (
            chemPointISAT<ThermoType>* phi0,
            const scalarField& phiq,
            scalarField& Rphiq
        );

        //- Check if the composition of the query point phiq lies in the
        //  ellipsoid of accuracy approximating the region of accuracy of the
        //  stored chemPoint phi0
        //  Input : phi0 the nearest chemPoint used in the linear interpolation
        //  phiq the composition of the query point for which we want to
        //  compute the mapping
        //  Output: true if phiq is in the EOA, false if not
        bool grow
        (
            chemPointISAT<ThermoType>* phi0,
            const scalarField& phiq,
            const scalarField& Rphiq
        );

        //- Clean and balance the tree
        bool cleanAndBalance();

        //- Functions to construct the gradients matrix
        //  When mechanism reduction is active, the A matrix is given by
        //        Aaa Aad
        //  A = ( Ada Add ), where the sub gradient matrices are:
        //  (Aaa) active species according to active species, (Aad) active
        //  species according to disabled species, (Ada) disabled species
        //  according to active species, and (Add) disabled species according to
        //  disabled species.
        //  The current implementation computes Aaa with the Jacobian of the
        //  reduced set of species. Aad = 0, Ada = 0, and Add = I.
        //  To be implemented: add options to compute the A matrix for different
        //  strategies
        void computeA
        (
            scalarSquareMatrix& A,
            const scalarField& Rphiq,
            const label li,
            const scalar rho,
            const scalar dt
        );


public:

    //- Runtime type information
    TypeName("ISAT");


    // Constructors

        //- Construct from dictionary
        ISAT
        (
            const dictionary& chemistryProperties,
            TDACChemistryModel<ThermoType>& chemistry
        );

        //- Disallow default bitwise copy construction
        ISAT(const ISAT&) = delete;


    // Destructor
    virtual ~ISAT();


    // Member Functions

        inline binaryTree<ThermoType>& chemisTree()
        {
            return chemisTree_;
        }

        inline const scalarField& scaleFactor() const
        {
            return scaleFactor_;
        }

        //- Return the size of the binary tree
        virtual inline label size()
        {
            return chemisTree_.size();
        }

        virtual void writePerformance();

        //- Find the closest stored leaf of phiQ and store the result in
        // RphiQ or return false.
        virtual bool retrieve
        (
            const Foam::scalarField& phiq,
            scalarField& Rphiq
        );

        //- Add information to the tabulation.
        //  This function can grow an existing point or add a new leaf to the
        //  binary tree Input : phiq the new composition to store Rphiq the
        //  mapping of the new composition point
        virtual label add
        (
            const scalarField& phiq,
            const scalarField& Rphiq,
            const label li,
            const scalar rho,
            const scalar deltaT
        );

        virtual bool update()
        {
            return cleanAndBalance();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryTabulationMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ISAT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

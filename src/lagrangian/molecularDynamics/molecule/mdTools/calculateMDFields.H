const List<DynamicList<molecule*>>& cellOccupancy = molecules.cellOccupancy();

forAll(cellOccupancy, cell)
{
    const List<molecule*>& molsInCell = cellOccupancy[cell];

    forAll(molsInCell, mIC)
    {
        molecule* mol = molsInCell[mIC];

        const label molId = mol->id();

        const vector& molU = mol->U();

        allSpeciesN_RU[molId][cell]++;

        allSpeciesM_RU[molId][cell] += mol->mass();

        allSpeciesVelocitySum_RU[molId][cell] += molU;

        allSpeciesVelocityMagSquaredSum_RU[molId][cell] += molU & molU;
    }
}

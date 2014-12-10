forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
{
    phaseModel& phase = iter();

    phase.DDtU() =
        fvc::ddt(phase.U())
      + fvc::div(phase.phi(), phase.U())
      - fvc::div(phase.phi())*phase.U();

    mrfZones.addCoriolis(phase.U(), phase.DDtU());
}

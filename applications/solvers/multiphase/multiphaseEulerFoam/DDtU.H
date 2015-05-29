forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
{
    phaseModel& phase = iter();

    phase.DDtU() =
        fvc::ddt(phase.U())
      + fvc::div(phase.phi(), phase.U())
      - fvc::div(phase.phi())*phase.U();

    MRF.addAcceleration(phase.U(), phase.DDtU());
}

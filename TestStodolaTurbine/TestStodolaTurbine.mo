
model TestStodolaTurbine

  WaterSteam.Machines.StodolaTurbine stodolaTurbine(Cst=2e6,eta_is_nom=0.94);
  WaterSteam.BoundaryConditions.SinkP puitsP(mode=0, P0=10000000);
  WaterSteam.BoundaryConditions.SourceP sourceP(option_temperature=2,mode=2,P0=27000000,h0=3475.e3);

equation

  connect(sourceP.C, stodolaTurbine.Ce);
  connect(stodolaTurbine.Cs, puitsP.C);

end TestStodolaTurbine;

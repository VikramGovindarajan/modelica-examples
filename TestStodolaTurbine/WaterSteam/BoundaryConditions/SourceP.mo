within WaterSteam.BoundaryConditions;

model SourceP "Water/steam source with fixed pressure"

  parameter Units.SI.AbsolutePressure P0=300000 "Source pressure";
  parameter Units.SI.Temperature T0=290 "Source temperature (active if option_temperature=1)";
  parameter Units.SI.SpecificEnthalpy h0=100000 "Source specific enthalpy (active if option_temperature=2)";
  parameter Integer option_temperature=1 "1:temperature fixed - 2:specific enthalpy fixed";
  parameter Integer mode=1 "IF97 region. 1:liquid - 2:steam - 4:saturation line - 0:automatic";

public
  Units.SI.AbsolutePressure P "Fluid pressure";
  Units.SI.MassFlowRate Q "Mass flow rate";
  Units.SI.Temperature T "Fluid temperature";
  Units.SI.SpecificEnthalpy h "Fluid enthalpy";
  ThermoSysPro.Properties.WaterSteam.Common.ThermoProperties_ph pro "Propriétés de l'eau";
  InstrumentationAndControl.Connectors.InputReal IPressure;
  InstrumentationAndControl.Connectors.InputReal ISpecificEnthalpy;
  Connectors.FluidOutlet C;
  InstrumentationAndControl.Connectors.InputReal ITemperature;

equation

  C.P = P;
  C.Q = Q;
  C.h_vol = h;

  if (cardinality(IPressure) == 0) then
    IPressure.signal = P0;
  end if;

  P = IPressure.signal;

  if (cardinality(ITemperature) == 0) then
      ITemperature.signal = T0;
  end if;

  if (cardinality(ISpecificEnthalpy) == 0) then
      ISpecificEnthalpy.signal = h0;
  end if;

  if (option_temperature == 1) then
    T = ITemperature.signal;
    h = ThermoSysPro.Properties.WaterSteam.IF97.SpecificEnthalpy_PT(P, T, 0);
  elseif (option_temperature == 2) then
    h = ISpecificEnthalpy.signal;
    T = pro.T;
  else
    assert(false, "SourcePressureWaterSteam: incorrect option");
  end if;

  pro = ThermoSysPro.Properties.WaterSteam.IF97.Water_Ph(P, h, mode);

end SourceP;

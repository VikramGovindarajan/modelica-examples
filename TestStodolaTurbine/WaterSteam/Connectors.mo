within WaterSteam;

package Connectors "Connectors"

  connector FluidInlet "Water/steam inlet fluid connector"
    Units.SI.AbsolutePressure P(start=1.e5)
      "Fluid pressure in the control volume";
    Units.SI.SpecificEnthalpy h_vol(start=1.e5)
      "Fluid specific enthalpy in the control volume";
    Units.SI.MassFlowRate Q(start=500)
      "Mass flow rate of the fluid crossing the boundary of the control volume";
    Units.SI.SpecificEnthalpy h(start=1.e5)
      "Specific enthalpy of the fluid crossing the boundary of the control volume";

    input Boolean a=true
      "Pseudo-variable for the verification of the connection orientation";
    output Boolean b
      "Pseudo-variable for the verification of the connection orientation";
  end FluidInlet;

  connector FluidInletI "Internal water/steam inlet fluid connector"
    Units.SI.AbsolutePressure P(start=1.e5)
      "Fluid pressure in the control volume";
    Units.SI.SpecificEnthalpy h_vol(start=1.e5)
      "Fluid specific enthalpy in the control volume";
    Units.SI.MassFlowRate Q(start=500)
      "Mass flow rate of the fluid crossing the boundary of the control volume";
    Units.SI.SpecificEnthalpy h(start=1.e5)
      "Specific enthalpy of the fluid crossing the boundary of the control volume";

    input Boolean a
      "Pseudo-variable for the verification of the connection orientation";
    output Boolean b
      "Pseudo-variable for the verification of the connection orientation";
  end FluidInletI;

  connector FluidOutlet "Water/steam outlet fluid connector"
    Units.SI.AbsolutePressure P(start=1.e5)
      "Fluid pressure in the control volume";
    Units.SI.SpecificEnthalpy h_vol(start=1.e5)
      "Fluid specific enthalpy in the control volume";
    Units.SI.MassFlowRate Q(start=500)
      "Mass flow rate of the fluid crossing the boundary of the control volume";
    Units.SI.SpecificEnthalpy h(start=1.e5)
      "Specific enthalpy of the fluid crossing the boundary of the control volume";

    output Boolean a
      "Pseudo-variable for the verification of the connection orientation";
    input Boolean b=true
      "Pseudo-variable for the verification of the connection orientation";

  end FluidOutlet;

  connector FluidOutletI "Internal water/steam outlet fluid connector"
    Units.SI.AbsolutePressure P(start=1.e5)
      "Fluid pressure in the control volume";
    Units.SI.SpecificEnthalpy h_vol(start=1.e5)
      "Fluid specific enthalpy in the control volume";
    Units.SI.MassFlowRate Q(start=500)
      "Mass flow rate of the fluid crossing the boundary of the control volume";
    Units.SI.SpecificEnthalpy h(start=1.e5)
      "Specific enthalpy of the fluid crossing the boundary of the control volume";

    output Boolean a
      "Pseudo-variable for the verification of the connection orientation";
    input Boolean b
      "Pseudo-variable for the verification of the connection orientation";

  end FluidOutletI;

end Connectors;

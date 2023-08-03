within Interfaces;

partial model PartialHeatTransfer "Common interface for heat transfer models"

  // Parameters
  replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
    "Medium in the component";
  
  import SI=Modelica.Units.SI;

  parameter Integer n=1 "Number of heat transfer segments";

  // Inputs provided to heat transfer model
  input Medium.ThermodynamicState[n] states
    "Thermodynamic states of flow segments";

  input SI.Area[n] surfaceAreas "Heat transfer areas";

  // Outputs defined by heat transfer model
  output SI.HeatFlowRate[n] Q_flows "Heat flow rates";

  // Parameters
  parameter Boolean use_k = false
    "= true to use k value for thermal isolation";
  parameter SI.CoefficientOfHeatTransfer k = 0
    "Heat transfer coefficient to ambient";
  parameter SI.Temperature T_ambient = system.T_ambient "Ambient temperature";
  outer System system "System wide properties";

  // Heat ports
  Interfaces.HeatPorts_a[n] heatPorts
    "Heat port to component boundary";

  // Variables
  SI.Temperature[n] Ts = Medium.temperature(states)
    "Temperatures defined by fluid states";

equation
  if use_k then
    Q_flows = heatPorts.Q_flow + {k*surfaceAreas[i]*(T_ambient - heatPorts[i].T) for i in 1:n};
  else
    Q_flows = heatPorts.Q_flow;
  end if;

end PartialHeatTransfer;

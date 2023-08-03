within Interfaces;

partial model PartialTwoPort "Partial component with two ports"

  import Modelica.Constants;

  outer System system "System wide properties";

  replaceable package Medium =
      Modelica.Media.Interfaces.PartialMedium "Medium in the component";

  parameter Boolean allowFlowReversal = system.allowFlowReversal
    "= true to allow flow reversal, false restricts to design direction (port_a -> port_b)";

  FluidPort_a port_a(
                                redeclare package Medium = Medium,
                     m_flow(min=if allowFlowReversal then -Constants.inf else 0))
    "Fluid connector a (positive design flow direction is from port_a to port_b)";

  FluidPort_b port_b(
                                redeclare package Medium = Medium,
                     m_flow(max=if allowFlowReversal then +Constants.inf else 0))
    "Fluid connector b (positive design flow direction is from port_a to port_b)";
  // Model structure, e.g., used for visualization

protected

  parameter Boolean port_a_exposesState = false
    "= true if port_a exposes the state of a fluid volume";

  parameter Boolean port_b_exposesState = false
    "= true if port_b.p exposes the state of a fluid volume";

  parameter Boolean showDesignFlowDirection = true
    "= false to hide the arrow in the model icon";

end PartialTwoPort;

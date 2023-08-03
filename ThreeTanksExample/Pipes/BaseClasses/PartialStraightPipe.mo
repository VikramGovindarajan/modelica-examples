within Pipes.BaseClasses;

partial model PartialStraightPipe "Base class for straight pipe models"
  extends Interfaces.PartialTwoPort;
  import SI=Modelica.Units.SI;
  
  // Geometry

  // Note: define nParallel as Real to support inverse calculations
  parameter Real nParallel(min=1)=1 "Number of identical parallel pipes";
  parameter SI.Length length "Length";
  parameter Boolean isCircular=true
    "= true, if cross sectional area is circular";
  parameter SI.Diameter diameter "Diameter of circular pipe";
  parameter SI.Area crossArea=Modelica.Constants.pi*diameter*diameter/4
    "Inner cross section area";
  parameter SI.Length perimeter(min=0)=Modelica.Constants.pi*diameter
    "Inner perimeter";
  parameter Types.Roughness roughness=2.5e-5
    "Average height of surface asperities (default: smooth steel pipe)";
  final parameter SI.Volume V=crossArea*length*nParallel "Volume size";

  // Static head
  parameter SI.Length height_ab=0 "Height(port_b) - Height(port_a)";

  // Pressure loss
  replaceable model FlowModel =
    Pipes.BaseClasses.FlowModels.DetailedPipeFlow
    constrainedby
    Pipes.BaseClasses.FlowModels.PartialStaggeredFlowModel
    "Wall friction, gravity, momentum flow";

equation
  assert(length >= height_ab, "Parameter length must be greater or equal height_ab.");

end PartialStraightPipe;

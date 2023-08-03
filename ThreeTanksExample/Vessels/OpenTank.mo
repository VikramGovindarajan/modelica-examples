within Vessels;

model OpenTank "Simple tank with inlet/outlet ports"
    import Modelica.Constants.pi;
	import SI=Modelica.Units.SI;
	import Types;

  // Tank properties
  SI.Height level(stateSelect=StateSelect.prefer, start=level_start_eps)
      "Level height of tank";
  SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";

  // Tank geometry
  parameter SI.Height height "Height of tank";
  parameter SI.Area crossArea "Area of tank";

  // Ambient
  parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure";
  parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature";

  // Initialization
  parameter SI.Height level_start(min=0) = 0.5*height
      "Start value of tank level";

  // Mass and energy balance, ports
  extends PartialLumpedVessel(
    final fluidVolume = V,
    final fluidLevel = level,
    final fluidLevel_max = height,
    final vesselArea = crossArea,
    heatTransfer(surfaceAreas={crossArea+2*sqrt(crossArea*pi)*level}),
    final initialize_p = false,
    final p_start = p_ambient);

  protected
  final parameter SI.Height level_start_eps = max(level_start, Modelica.Constants.eps);

equation
  // Total quantities
  V = crossArea*level "Volume of fluid";
  medium.p = p_ambient;

  // Source termsEnergy balance
  if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
    Wb_flow = 0
        "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
  else
    Wb_flow = -p_ambient*der(V);
  end if;

  //Determine port properties
  for i in 1:nPorts loop
    vessel_ps_static[i] = max(0, level - portsData_height[i])*system.g*medium.d + p_ambient;
  end for;

initial equation
  if massDynamics == Types.Dynamics.FixedInitial then
    level = level_start_eps;
  elseif massDynamics == Types.Dynamics.SteadyStateInitial then
    der(level) = 0;
  end if;

end OpenTank;

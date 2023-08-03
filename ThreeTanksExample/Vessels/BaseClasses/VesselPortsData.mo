within Vessels.BaseClasses;

record VesselPortsData "Data to describe inlet/outlet ports at vessels:
diameter -- Inner (hydraulic) diameter of inlet/outlet port
height -- Height over the bottom of the vessel
zeta_out -- Hydraulic resistance out of vessel, default 0.5 for small diameter mounted flush with the wall
zeta_in -- Hydraulic resistance into vessel, default 1.04 for small diameter mounted flush with the wall"
      // extends Modelica.Icons.Record;
  import SI=Modelica.Units.SI;
  
  parameter SI.Diameter diameter
    "Inner (hydraulic) diameter of inlet/outlet port";
  parameter SI.Height height = 0 "Height over the bottom of the vessel";
  parameter Real zeta_out(min=0)=0.5
    "Hydraulic resistance out of vessel, default 0.5 for small diameter mounted flush with the wall";
  parameter Real zeta_in(min=0)=1.04
    "Hydraulic resistance into vessel, default 1.04 for small diameter mounted flush with the wall";
  
end VesselPortsData;

within InstrumentationAndControl;

package Connectors "Connectors"


  connector InputReal
    input Real signal;

  end InputReal;

  connector InputLogical
     input Boolean signal;
     annotation (
       Window(
         x=0.37,
         y=0.02,
         width=0.49,
         height=0.65),
       Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{-100,100},{-100,-100},{100,0},{-100,100}},
            lineColor={0,0,255},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid)}),
       Documentation(info="<html>
<p><b>Version 1.0</b></p>
</HTML>
"));
  end InputLogical;

  connector InputInteger
    input Integer signal;
    annotation (
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics),
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{-100,100},{-100,-100},{100,0},{-100,100}},
            lineColor={0,0,255},
            fillColor={255,0,255},
            fillPattern=FillPattern.Solid)}),
      Window(
        x=0.34,
        y=0.2,
        width=0.6,
        height=0.6),
      Documentation(info="<html>
<p><b>Version 1.6</b></p>
</HTML>
"));

  end InputInteger;

  connector OutputInteger
    output Integer signal;
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{-100,100},{-100,-100},{100,0},{-100,100}},
            lineColor={0,0,255},
            fillColor={255,0,128},
            fillPattern=FillPattern.Solid)}),
      Window(
        x=0.34,
        y=0.18,
        width=0.6,
        height=0.6),
      Documentation(info="<html>
<p><b>Version 1.6</b></p>
</HTML>
"));
  end OutputInteger;

  connector OutputLogical
      output Boolean signal;
      annotation (
        Window(
          x=0.29,
          y=0.11,
          width=0.6,
          height=0.6),
        Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{-100,100},{-100,-100},{100,0},{-100,100}},
            lineColor={0,0,255},
            fillColor={127,255,0},
            fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
<p><b>Version 1.0</b></p>
</HTML>
"));
  end OutputLogical;

  connector OutputReal
      output Real signal;
      annotation (
        Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{-100,100},{-100,-100},{100,0},{-100,100}},
            lineColor={0,0,255},
            fillColor={0,255,255},
            fillPattern=FillPattern.Solid)}),
        Window(
          x=0.34,
          y=0.18,
          width=0.6,
          height=0.6),
        Documentation(info="<html>
<p><b>Version 1.0</b></p>
</HTML>
"));
  end OutputReal;
  annotation (
    Window(
      x=0.05,
      y=0.26,
      width=0.25,
      height=0.25,
      library=1,
      autolayout=1),
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2},
        initialScale=0.1), graphics={
        Text(
          extent={{-102,0},{24,-26}},
          lineColor={242,148,0},
          textString=
               "Thermo"),
        Text(
          extent={{-4,8},{68,-34}},
          lineColor={46,170,220},
          textString=
               "SysPro"),
        Polygon(
          points={{-62,2},{-58,4},{-48,8},{-32,12},{-16,14},{6,14},{26,12},{42,
              8},{52,2},{42,6},{28,10},{6,12},{-12,12},{-16,12},{-34,10},{-50,6},
              {-62,2}},
          lineColor={46,170,220},
          fillColor={46,170,220},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-44,38},{-24,38},{-26,30},{-26,22},{-24,14},{-24,12},{-46,8},
              {-42,22},{-42,30},{-44,38}},
          lineColor={46,170,220},
          fillColor={46,170,220},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-26,20},{-20,20},{-20,22},{-14,22},{-14,20},{-12,20},{-12,12},
              {-26,12},{-28,12},{-26,20}},
          lineColor={46,170,220},
          fillColor={46,170,220},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-8,14},{-8,24},{-6,24},{-6,14},{-8,14}},
          lineColor={46,170,220},
          fillColor={46,170,220},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-8,30},{-6,26}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-8,36},{-6,32}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-8,42},{-6,38}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-8,48},{-6,44}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-4,14},{-4,26},{-2,26},{-2,14},{-4,14}},
          lineColor={46,170,220},
          fillColor={46,170,220},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-4,32},{-2,28}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-4,38},{-2,34}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-4,44},{-2,40}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-4,50},{-2,46}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-2,20},{8,20},{8,22},{10,22},{18,22},{18,12},{-4,14},{-2,20}},
          lineColor={46,170,220},
          fillColor={46,170,220},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-62,2},{-58,4},{-48,8},{-36,10},{-18,12},{6,12},{26,10},{42,
              6},{52,0},{42,4},{28,8},{6,10},{-12,10},{-18,10},{-38,8},{-50,6},
              {-62,2}},
          lineColor={242,148,0},
          fillColor={242,148,0},
          fillPattern=FillPattern.Solid),
        Line(
          points={{22,12},{22,14},{22,16},{24,14},{20,18}},
          color={46,170,220},
          thickness=0.5),
        Line(
          points={{26,12},{26,14},{26,16},{28,14},{24,18}},
          color={46,170,220},
          thickness=0.5),
        Line(
          points={{30,10},{30,12},{30,14},{32,12},{28,16}},
          color={46,170,220},
          thickness=0.5),
        Polygon(
          points={{36,8},{36,30},{34,34},{36,38},{40,38},{40,8},{36,8}},
          lineColor={46,170,220},
          fillColor={46,170,220},
          fillPattern=FillPattern.Solid),
        Rectangle(extent={{-100,80},{80,-100}}, lineColor={0,0,255}),
        Line(
          points={{-100,80},{-80,100},{100,100},{100,-80},{80,-100}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{80,80},{100,100}},
          color={0,0,255},
          smooth=Smooth.None)}),
    Documentation(info="<html>
<p><b>Version 1.0</b></p>
</HTML>
"));

end Connectors;
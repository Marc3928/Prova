within Flow1D2PHD_Library.FluidPH;

model Storage
    replaceable package Medium = Modelica.Media.Water.WaterIF97_ph constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model";
    replaceable model HeatTransfer = Flow1D2PHD_Library.HeatTransfer.Constant_gamma (redeclare package Medium=Medium, L=L , Dhyd=Dhyd, gamma = 100) constrainedby Flow1D2PHD_Library.HeatTransfer.BaseClasses.BaseHeatTransfer; 
    Interfaces.Flange_wlwv_A port_A (redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {-98, 0}, extent = {{-10, -10}, {10, 10}})));
    Interfaces.Flange_wlwv_B port_B (redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {98, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {98, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heat_port annotation(
    Placement(transformation(origin = {0, 96}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {0, 96}, extent = {{-10, -10}, {10, 10}})));
    Medium.ThermodynamicState fluidState;
    
    
    HeatTransfer heatTransfer(w=wht, fluidState=fluidState);  
    parameter Modelica.Units.SI.Length Dhyd;
    parameter Modelica.Units.SI.Length L;
    parameter Modelica.Units.SI.Volume V;
    Modelica.Units.SI.MassFlowRate wht;
    
    Medium.ThermodynamicState state_in;
    Medium.SaturationProperties sat_in;
    Medium.Density rhov_in "Densità di saturazione vapore";
    Medium.Density rhol_in "Densità di saturazione liquido";
    Medium.Density rho_in "Densità media della miscela bifase";
    Medium.SpecificEnthalpy hv_in "Entalpia di saturazione vapore";
    Medium.SpecificEnthalpy hl_in "Entalpia di saturazione liquido";
    Modelica.Units.SI.QualityFactor x_in "Titolo di vapore";
    Modelica.Units.SI.QualityFactor gamma_in "Grado di vuoto";
    Modelica.Units.SI.QualityFactor flow_quality_in "Flow quality";
    Modelica.Units.SI.MassFlowRate w_in "Portata massica";
    
    Modelica.Units.SI.MassFlowRate w_L_A "Portata massica vapore alle interfacce";
    Modelica.Units.SI.MassFlowRate w_V_A "Portata massica vapore alle interfacce";
    Modelica.Units.SI.MassFlowRate w_L_B "Portata massica liquido alle interfacce";
    Modelica.Units.SI.MassFlowRate w_V_B "Portata massica vapore alle interfacce";
    
    Modelica.Units.SI.SpecificEnthalpy h_L_A "Portata massica vapore alle interfacce";
    Modelica.Units.SI.SpecificEnthalpy h_V_A "Portata massica vapore alle interfacce";
    Modelica.Units.SI.SpecificEnthalpy h_L_B "Portata massica liquido alle interfacce";
    Modelica.Units.SI.SpecificEnthalpy h_V_B "Portata massica vapore alle interfacce";  
  
    Medium.Density rho "Densità media della miscela bifase";
    Medium.Density rhov "Densità di saturazione vapore";
    Medium.Density rhol "Densità di saturazione liquido"; 
    Medium.ThermodynamicState state;
    Medium.SaturationProperties sat;
    Medium.Temperature Ts "Temperatura di saturazione";
    Medium.SpecificEnthalpy hv "Entalpia di saturazione vapore";
    Medium.SpecificEnthalpy hl "Entalpia di saturazione liquido";
    Medium.DerDensityByEnthalpy drdh; 
    Medium.DerDensityByPressure drdp;
    Modelica.Units.SI.QualityFactor x "Titolo di vapore";
    Modelica.Units.SI.QualityFactor alpha "Grado di vuoto";
    Modelica.Units.SI.Mass M_L "Massa di liquido";
    Modelica.Units.SI.Mass M_V "Massa di vapore";
    
    
    // Variabili di stato
    Real p( start = 900000.0,  fixed = true) "Pressione nel volume j (Pa)"; 
    Real h( start = 1100000.0,  fixed = true) "Entalpia specifica nel volume j (J/kg)";
    
    //Definire pstart - hstart
     
equation
    connect (heat_port, heatTransfer.port);
    
    state = Medium.setState_ph(p, h);
    sat = Medium.setSat_p(p);
    Ts   = Medium.saturationTemperature(p);
    rhov = Medium.dewDensity(sat);
    rhol = Medium.bubbleDensity(sat);
    hv = Medium.dewEnthalpy(sat);
    hl = Medium.bubbleEnthalpy(sat);
    drdh = Medium.density_derh_p(state);
    drdp = Medium.density_derp_h(state);
    x = Medium.vapourQuality(state);
    rho = 1/((x/rhov)+((1-x)/rhol));
    alpha = (rho-rhol)/(rhov-rhol);
    M_L = (1-x)*rho*V;
    M_V = x*rho*V;
    wht = 0.25 * (abs(w_L_A) + abs(w_L_B) + abs(w_V_A) + abs(w_V_B));
    
    // Massa
    V * (drdp * der(p) + drdh * der(h)) = (w_L_A + w_V_A) + (w_L_B + w_V_B);
    // Energia
    V*rho*der(h) = (w_V_A * h_V_A + w_L_A * h_L_A) + (w_V_B * h_V_B + w_L_B * h_L_B) + V*der(p) + heatTransfer.Q_flow;
    
    w_L_A = port_A.liq.m_flow;
    w_V_A = port_A.vap.m_flow;
    w_L_B = port_B.liq.m_flow;
    w_V_B = port_B.vap.m_flow;
    
    h_L_A = actualStream(port_A.liq.h_outflow);
    h_V_A = actualStream(port_A.vap.h_outflow);
    h_L_B = actualStream(port_B.liq.h_outflow);
    h_V_B = actualStream(port_B.vap.h_outflow);
    
    port_A.liq.p = p;
    port_A.vap.p = p;
    port_B.liq.p = p;
    port_B.vap.p = p;
    
    port_A.liq.h_outflow = hl;
    port_A.vap.h_outflow = hv;
    port_B.liq.h_outflow = hl;
    port_B.vap.h_outflow = hv;
    
    port_A.liq.Xi_outflow = inStream(port_B.liq.Xi_outflow);
    port_A.vap.Xi_outflow = inStream(port_B.vap.Xi_outflow);
    port_B.liq.Xi_outflow = inStream(port_A.liq.Xi_outflow);
    port_B.vap.Xi_outflow = inStream(port_A.vap.Xi_outflow);
    
    port_A.liq.C_outflow = inStream(port_B.liq.C_outflow);
    port_B.liq.C_outflow = inStream(port_A.liq.C_outflow);
    port_A.vap.C_outflow = inStream(port_B.vap.C_outflow);
    port_B.vap.C_outflow = inStream(port_A.vap.C_outflow);

annotation(
    Icon(graphics = {Ellipse(fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
end Storage;

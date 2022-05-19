/* Tank emptying code extracts */
/* (c) AspireSpace */
/* Variables are in metric units except where stated otherwise */
#include "nitrous_oxide.h" /* header file for nitrous oxide property calcs subroutines */
/* prototypes */
static double injector_model(double upstream_pressure, double downstream_pressure);
/* square */
double SQR(double bob)
{
	return((bob * bob));
}
#define CENTIGRADE_TO_KELVIN 273.15 // to Kelvin
#define BAR_TO_PASCALS 100000.0
#define PASCALS_TO_BAR (1.0 / BAR_TO_PASCALS)
/* calculate injector pressure drop (Bar) and mass flowrate (kg/sec) */
static double injector_model(double upstream_pressure, double downstream_pressure)
{
	double mass_flowrate;
	double pressure_drop;
	pressure_drop = upstream_pressure - downstream_pressure; /* Bar */
															 /* reality check */
	if (pressure_drop < 0.00001)
		pressure_drop = 0.00001;
	/* is injector pressure drop lower than 20 percent of chamber pressure? */
	if ((pressure_drop / hybrid.chamber_pressure_bar) < 0.2)
		hybrid.hybrid_fault = 3; // too low for safety
								 /* Calculate fluid flowrate through the injector, based on the */
								 /* total-pressure loss factor between the tank and combustion chamber */
								 /* (injector_loss_coefficient includes K coefficient and orifice cross-sectional areas) */
	mass_flowrate = sqrt((2.0 * hybrid.tank_liquid_density * pressure_drop / hybrid.injector_loss_coefficient));

	return(mass_flowrate); /* kg/sec */
}
/* Equilibrium (instantaneous boiling) tank blowdown model */
void Nitrous_tank(void)
{
	double bob;
	double Chamber_press_bar_abs;
	double delta_outflow_mass, deltaQ, deltaTemp;
	double Enth_of_vap;
	double Spec_heat_cap;
	double tc;
	static double lagged_bob = 0.0;
	/* blowdown simulation using nitrous oxide property calcs subroutines */

	/* update last-times values, O = 'old' */
	Omdot_tank_outflow = mdot_tank_outflow;
	Enth_of_vap = nox_enthV(hybrid.tank_fluid_temperature_K); /* Get enthalpy (latent heat) of vaporisation */
	Spec_heat_cap = nox_CpL(hybrid.tank_fluid_temperature_K); /* Get specific heat capacity of the liquid nitrous */
															  /* Calculate the heat removed from the liquid nitrous during its vaporisation */
	deltaQ = vaporised_mass_old * Enth_of_vap;
	/* temperature drop of the remaining liquid nitrous due to losing this heat */
	deltaTemp = -(deltaQ / (hybrid.tank_liquid_mass * Spec_heat_cap));
	hybrid.tank_fluid_temperature_K += deltaTemp; /* update fluid temperature */
												  /* reality checks */
	if (hybrid.tank_fluid_temperature_K < (-90.0 + CENTIGRADE_TO_KELVIN))
	{
		hybrid.tank_fluid_temperature_K = (-90.0 + CENTIGRADE_TO_KELVIN); /* lower limit */
		hybrid.hybrid_fault = 1;
	}
	else if (hybrid.tank_fluid_temperature_K >(36.0 + CENTIGRADE_TO_KELVIN))
	{
		hybrid.tank_fluid_temperature_K = (36.0 + CENTIGRADE_TO_KELVIN); /* upper limit */
		hybrid.hybrid_fault = 2;
	}
	/* get current nitrous properties */
	hybrid.tank_liquid_density = nox_Lrho(hybrid.tank_fluid_temperature_K);
	hybrid.tank_vapour_density = nox_Vrho(hybrid.tank_fluid_temperature_K);
	hybrid.tank_pressure_bar = nox_vp(hybrid.tank_fluid_temperature_K); /* vapour pressure, Bar abs */
	Chamber_press_bar_abs = hybrid.chamber_pressure_bar; /* Bar Abs */

														 /* calculate injector pressure drop and mass flowrate */
	mdot_tank_outflow = injector_model(hybrid.tank_pressure_bar, Chamber_press_bar_abs);
	/* integrate mass flowrate using Addams second order integration formula */
	/* (my preferred integration formulae, feel free to choose your own.) */
	/* Xn=X(n-1) + DT/2 * ((3 * Xdot(n-1) - Xdot(n-2)) */
	/* O infront of a variable name means value from previous timestep (Old) */
	delta_outflow_mass = 0.5 * delta_time * (3.0 * mdot_tank_outflow - Omdot_tank_outflow);
	/* drain the tank based on flowrates only */
	hybrid.tank_propellant_contents_mass -= delta_outflow_mass; /* update mass within tank for next iteration */
	old_liquid_nox_mass -= delta_outflow_mass; /* update liquid mass within tank for next iteration */

											   /* now the additional effects of phase changes */
											   /* The following equation is applicable to the nitrous tank, containing saturated nitrous: */
											   /* tank_volume = liquid_nox_mass / liquid_nox_density + gaseous_nox_mass / gaseous_nox_density */
											   /* Rearrage this equation to calculate current liquid_nox_mass */
	bob = (1.0 / hybrid.tank_liquid_density) - (1.0 / hybrid.tank_vapour_density);
	hybrid.tank_liquid_mass
		= (hybrid.tank_volume - (hybrid.tank_propellant_contents_mass / hybrid.tank_vapour_density)) / bob;
	hybrid.tank_vapour_mass = hybrid.tank_propellant_contents_mass - hybrid.tank_liquid_mass;
	/* update for next iteration */
	bob = old_liquid_nox_mass - hybrid.tank_liquid_mass;

	/* Add a 1st-order time lag (of 0.15 seconds) to aid numerical */
	/* stability (this models the finite time required for boiling) */
	tc = delta_time / 0.15;
	lagged_bob = tc * (bob - lagged_bob) + lagged_bob; // 1st-order lag
	vaporised_mass_old = lagged_bob;
	/* update tank contents for next iteration */
	old_liquid_nox_mass = hybrid.tank_liquid_mass;
}
/* Subroutine to initialise program variables */
/* Gets called only once, prior to firing */
void initialise_hybrid_engine(void)
{
	hybrid.hybrid_fault = 0;
	hybrid.tank_vapour_mass = 0.0;
	mdot_tank_outflow = 0.0;

	/* set either initial nitrous vapour (tank) pressure */
	/* or initial nitrous temperature (deg Kelvin) */
	if (hybrid.press_or_temp == true)
		hybrid.tank_fluid_temperature_K = nox_on_press(hybrid.initial_tank_pressure); // set tank pressure
	else
	{ /* set nitrous temperature */
		hybrid.tank_fluid_temperature_K = hybrid.initial_fluid_propellant_temp + CENTIGRADE_TO_KELVIN;
		hybrid.initial_tank_pressure = nox_vp(hybrid.tank_fluid_temperature_K);
	}
	/* reality check */
	if (hybrid.tank_fluid_temperature_K > (36.0 + CENTIGRADE_TO_KELVIN))
	{
		hybrid.tank_fluid_temperature_K = 36.0 + CENTIGRADE_TO_KELVIN;
		hybrid.hybrid_fault = 2;
	}
	/* get initial nitrous properties */
	hybrid.tank_liquid_density = nox_Lrho(hybrid.tank_fluid_temperature_K);
	hybrid.tank_vapour_density = nox_Vrho(hybrid.tank_fluid_temperature_K);
	/* base the nitrous vapour volume on the tank percentage ullage (gas head-space) */
	hybrid.tank_vapour_volume = (hybrid.initial_ullage / 100.0) * hybrid.tank_volume;
	hybrid.tank_liquid_volume = hybrid.tank_volume - hybrid.tank_vapour_volume;
	hybrid.tank_liquid_mass = hybrid.tank_liquid_density * hybrid.tank_liquid_volume;
	hybrid.tank_vapour_mass = hybrid.tank_vapour_density * hybrid.tank_vapour_volume;
	hybrid.tank_propellant_contents_mass = hybrid.tank_liquid_mass + hybrid.tank_vapour_mass; /* total mass within tank */
																							  /* initialise values needed later */
	old_liquid_nox_mass = hybrid.tank_liquid_mass;
	old_vapour_nox_mass = hybrid.tank_vapour_mass;
	hybrid.initial_liquid_propellant_mass = hybrid.tank_liquid_mass;
	hybrid.initial_vapour_propellant_mass = hybrid.tank_vapour_mass;
	/* guessed initial value of amount of nitrous vaporised per iteration */
	/* in the nitrous tank blowdown model (actual value is not important) */
	vaporised_mass_old = 0.001;

	/* individual injector orifice total loss coefficent K2 */
	bob = PI * SQR((hybrid.orifice_diameter / 2.0)); /* orifice cross sectional area */
	hybrid.injector_loss_coefficient
		= (hybrid.orifice_k2_coefficient / (SQR((hybrid.orifice_number * bob)))) * PASCALS_TO_BAR;
}
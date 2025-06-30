
////////////////////////////////////////////////////////////////////////////
//
// accrete planet system generator
//
// extended version of accrete "C" source code
//
// modified 2025.06.30 version number used here 0000.0000.0002
//
// based on
//
// "C version, by Matt Burdick, 1988. Most commonly distributed as accrete."
//
// https://znark.com/create/files/accrete.zip
//
// https://znark.com/create/accrete.html
//
//////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------*/
/*                           BIBLIOGRAPHY                               */
/*  Dole, Stephen H.  "Formation of Planetary Systems by Aggregation:   */
/*      a Computer Simulation"  October 1969,  Rand Corporation Paper   */
/*	P-4226.								*/
/*----------------------------------------------------------------------*/

//    Dole, Stephen H. "Computer Simulation of the Formation of Planetary Systems". Icarus, vol 13, pp 494-508, 1970.
//    Isaacman, Richard. & Sagan, Carl. "Computer Simulation of Planetary Accretion Dynamics: Sensitivity to Initial Conditions". Icarus, vol 31, p 510, 1977.
//    Fogg, Martyn J. "Extra-Solar Planetary Systems: A Microcomputer Simulation". Journal of the British // Interplanetary Society, vol 38, p 501-514, 1985.


//
// compilation with "accreteb.h" on linux gcc 
//
// gcc accreteb.c -lm -o accreteb
//


#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // For strcmp (string comparison)
#include <time.h>
#include <math.h>
#include "accreteb.h"



double radians_per_rotation = 6.283185307179586;

int seed1=777;

int render=1;


double stellar_mass_ratio = 1;
double stellar_luminosity_ratio = 1.0;
double main_seq_life = 1.0e10;
double age = 4.5e9;
double r_ecosphere = 1.0;
int spin_resonance = 0;

double base_dust_limit=200;
double base_innermost_planet=0.3;
double base_stellar_mass_ratio=50.0; ## outermost planet is prop sqrt of this
double base_cloud_ecentricity=0.2;
// simple migration: change distance of planet
int use_migration = 0;  
int use_coeff_migration =1;
int use_exp_migration = 0;   
int use_random_migration =0;
double migration_coeff =0.1; // for migration bu coeff
double migration_exp =2.5; // for migration by randim



// Prototyypit puuttuville funktioille
double random_number(double inner, double outer);
double random_eccentricity(void);
double about(double value, double variation);
void apply_migration(planet_pointer planet);
int generate_and_render_povray(void);
int display_system(void);

typedef struct dust_bands_record  *dust_pointer;
typedef struct dust_bands_record {
	  double inner_edge;
	  double outer_edge;
	  int dust_present;
	  int gas_present;
	  dust_pointer next_band;
     } dust_bands;

/* A few variables global to the entire program:		*/
extern planet_pointer planet_head;

/* Now for some variables global to the accretion process:      */
int dust_left;
double r_inner, r_outer, reduced_mass, dust_density, cloud_eccentricity;
dust_pointer dust_head;


void set_initial_conditions(inner_limit_of_dust, outer_limit_of_dust)
double inner_limit_of_dust, outer_limit_of_dust;
{
     dust_head = (dust_bands *)malloc(sizeof(dust_bands));
     planet_head = NULL;
     dust_head->next_band = NULL;
     dust_head->outer_edge = outer_limit_of_dust;
     dust_head->inner_edge = inner_limit_of_dust;
     dust_head->dust_present = TRUE;
     dust_head->gas_present = TRUE;
     dust_left = TRUE;
     cloud_eccentricity = base_cloud_ecentricity;
}




double stellar_dust_limit(stellar_mass_ratio)
double stellar_mass_ratio;
{
     return(base_dust_limit * pow(stellar_mass_ratio,(1.0 / 3.0)));
}

double innermost_planet(stellar_mass_ratio)
double stellar_mass_ratio;
{
     return(base_innermost_planet * pow(stellar_mass_ratio,(1.0 / 3.0)));
}

double outermost_planet(stellar_mass_ratio)
double stellar_mass_ratio;
{
     return(base_stellar_mass_ratio * pow(stellar_mass_ratio,(1.0 / 3.0)));
}

double inner_effect_limit(a, e, mass)
double a, e, mass;
{
     return (a * (1.0 - e) * (1.0 - mass) / (1.0 + cloud_eccentricity));
}

double outer_effect_limit(a, e, mass)
double a, e, mass;
{
     return (a * (1.0 + e) * (1.0 + reduced_mass) / (1.0 - cloud_eccentricity));
}

int dust_available(inside_range, outside_range)
double inside_range, outside_range;
{
     dust_pointer current_dust_band;
     int dust_here;

     current_dust_band = dust_head;
     while ((current_dust_band != NULL)
	    && (current_dust_band->outer_edge < inside_range))
	  current_dust_band = current_dust_band->next_band;
     if (current_dust_band == NULL)
	  dust_here = FALSE;
     else dust_here = current_dust_band->dust_present;
     while ((current_dust_band != NULL)
	    && (current_dust_band->inner_edge < outside_range)) {
	       dust_here = dust_here || current_dust_band->dust_present;
	       current_dust_band = current_dust_band->next_band;
	  }
     return(dust_here);
}

void update_dust_lanes(min, max, mass, crit_mass,
		       body_inner_bound, body_outer_bound)
double min, max, mass, crit_mass, body_inner_bound, body_outer_bound;
{
     int gas;
     dust_pointer node1, node2, node3;

     dust_left = FALSE;
     if ((mass > crit_mass))
	  gas = FALSE;
     else
	  gas = TRUE;
     node1 = dust_head;
     while ((node1 != NULL))
     {
	  if (((node1->inner_edge < min) && (node1->outer_edge > max)))
	  {
	       node2 = (dust_bands *)malloc(sizeof(dust_bands));
	       node2->inner_edge = min;
	       node2->outer_edge = max;
	       if ((node1->gas_present == TRUE))
		    node2->gas_present = gas;
	       else
		    node2->gas_present = FALSE;
	       node2->dust_present = FALSE;
	       node3 = (dust_bands *)malloc(sizeof(dust_bands));
	       node3->inner_edge = max;
	       node3->outer_edge = node1->outer_edge;
	       node3->gas_present = node1->gas_present;
	       node3->dust_present = node1->dust_present;
	       node3->next_band = node1->next_band;
	       node1->next_band = node2;
	       node2->next_band = node3;
	       node1->outer_edge = min;
	       node1 = node3->next_band;
	  }
	  else
	       if (((node1->inner_edge < max) && (node1->outer_edge > max)))
	       {
		    node2 = (dust_bands *)malloc(sizeof(dust_bands));
		    node2->next_band = node1->next_band;
		    node2->dust_present = node1->dust_present;
		    node2->gas_present = node1->gas_present;
		    node2->outer_edge = node1->outer_edge;
		    node2->inner_edge = max;
		    node1->next_band = node2;
		    node1->outer_edge = max;
		    if ((node1->gas_present == TRUE))
			 node1->gas_present = gas;
		    else
			 node1->gas_present = FALSE;
		    node1->dust_present = FALSE;
		    node1 = node2->next_band;
	       }
	       else
		    if (((node1->inner_edge < min) && (node1->outer_edge > min)))
		    {
			 node2 = (dust_bands *)malloc(sizeof(dust_bands));
			 node2->next_band = node1->next_band;
			 node2->dust_present = FALSE;
			 if ((node1->gas_present == TRUE))
			      node2->gas_present = gas;
			 else
			      node2->gas_present = FALSE;
			 node2->outer_edge = node1->outer_edge;
			 node2->inner_edge = min;
			 node1->next_band = node2;
			 node1->outer_edge = min;
			 node1 = node2->next_band;
		    }
		    else
			 if (((node1->inner_edge >= min) && (node1->outer_edge <= max)))
			 {
			      if ((node1->gas_present == TRUE))
				   node1->gas_present = gas;
			      node1->dust_present = FALSE;
			      node1 = node1->next_band;
			 }
			 else
			      if (((node1->outer_edge < min) || (node1->inner_edge > max)))
				   node1 = node1->next_band;
     }
     node1 = dust_head;
     while ((node1 != NULL))
     {
	  if (((node1->dust_present)
	       && (((node1->outer_edge >= body_inner_bound)
		    && (node1->inner_edge <= body_outer_bound)))))
	       dust_left = TRUE;
	  node2 = node1->next_band;
	  if ((node2 != NULL))
	  {
	       if (((node1->dust_present == node2->dust_present)
		    && (node1->gas_present == node2->gas_present)))
	       {
		    node1->outer_edge = node2->outer_edge;
		    node1->next_band = node2->next_band;
		    free(node2);
	       }
	  }
	  node1 = node1->next_band;
     }
}

double collect_dust(last_mass, a, e, crit_mass, dust_band)
double last_mass, a, e, crit_mass;
dust_pointer dust_band;
{
     double mass_density, temp1, temp2, temp, temp_density, bandwidth, width, volume;

     temp = last_mass / (1.0 + last_mass);
     reduced_mass = pow(temp,(1.0 / 4.0));
     r_inner = inner_effect_limit(a, e, reduced_mass);
     r_outer = outer_effect_limit(a, e, reduced_mass);
     if ((r_inner < 0.0))
	  r_inner = 0.0;
     if ((dust_band == NULL))
	  return(0.0);
     else
     {
	  if ((dust_band->dust_present == FALSE))
	       temp_density = 0.0;
	  else
	       temp_density = dust_density;
	  if (((last_mass < crit_mass) || (dust_band->gas_present == FALSE)))
	       mass_density = temp_density;
	  else
	       mass_density = K * temp_density / (1.0 + sqrt(crit_mass / last_mass)
						  * (K - 1.0));
	  if (((dust_band->outer_edge <= r_inner)
	       || (dust_band->inner_edge >= r_outer)))
	       return(collect_dust(last_mass,a,e,crit_mass, dust_band->next_band));
	  else
	  {
	       bandwidth = (r_outer - r_inner);
	       temp1 = r_outer - dust_band->outer_edge;
	       if (temp1 < 0.0)
		    temp1 = 0.0;
	       width = bandwidth - temp1;
	       temp2 = dust_band->inner_edge - r_inner;
	       if (temp2 < 0.0)
		    temp2 = 0.0;
	       width = width - temp2;
	       temp = 4.0 * PI * pow(a,2.0) * reduced_mass
		    * (1.0 - e * (temp1 - temp2) / bandwidth);
	       volume = temp * width;
	       return(volume * mass_density
		      + collect_dust(last_mass,a,e,crit_mass,
				     dust_band->next_band));
	  }
     }
}


/*--------------------------------------------------------------------------*/
/*   Orbital radius is in AU, eccentricity is unitless, and the stellar     */
/*  luminosity ratio is with respect to the sun.  The value returned is the */
/*  mass at which the planet begins to accrete gas as well as dust, and is  */
/*  in units of solar masses.                                               */
/*--------------------------------------------------------------------------*/

double critical_limit(orbital_radius, eccentricity, stellar_luminosity_ratio)
double orbital_radius, eccentricity, stellar_luminosity_ratio;
{
     double temp, perihelion_dist;

     perihelion_dist = (orbital_radius - orbital_radius * eccentricity);
     temp = perihelion_dist * sqrt(stellar_luminosity_ratio);
     return(B * pow(temp,-0.75));
}



void accrete_dust(seed_mass, a, e, crit_mass,
		  body_inner_bound, body_outer_bound)
double *seed_mass, a, e, crit_mass,
     body_inner_bound, body_outer_bound;
{
     double perihelion_dist, new_mass, temp_mass;

     new_mass = (*seed_mass);
     do
     {
	  temp_mass = new_mass;
	  new_mass = collect_dust(new_mass,a,e,crit_mass,
				  dust_head);
     }
     while (!(((new_mass - temp_mass) < (0.0001 * temp_mass))));
     (*seed_mass) = (*seed_mass) + new_mass;
     update_dust_lanes(r_inner,r_outer,(*seed_mass),crit_mass,body_inner_bound,body_outer_bound);
}



void coalesce_planetesimals(double a, double e, double mass, double crit_mass,
                             double stellar_luminosity_ratio,
                             double body_inner_bound, double body_outer_bound)
{
    planet_pointer current = planet_head;
    planet_pointer previous = NULL;
    int merged = FALSE;

    while (current != NULL)
    {
        double delta_a = current->a - a;
        double reduced_mass = pow(current->mass / (1.0 + current->mass), 0.25);
        double reach_a = (delta_a > 0.0)
                         ? (a * (1.0 + e) * (1.0 + reduced_mass)) - a
                         : a - (a * (1.0 - e) * (1.0 - reduced_mass));
        double reach_current = (delta_a > 0.0)
                               ? current->a - (current->a * (1.0 - current->e) * (1.0 - reduced_mass))
                               : (current->a * (1.0 + current->e) * (1.0 + reduced_mass)) - current->a;

        if (fabs(delta_a) <= fmax(fabs(reach_a), fabs(reach_current)))
        {
#ifdef VERBOSE
            printf("Collision between planetesimals at %.3lf AU and %.3lf AU\n", a, current->a);
#endif
            double new_a = (current->mass + mass) /
                           ((current->mass / current->a) + (mass / a));

            double total_angular_momentum = 
                current->mass * sqrt(current->a) * sqrt(1.0 - pow(current->e, 2.0)) +
                mass * sqrt(a) * sqrt(1.0 - pow(e, 2.0));

            double new_ecc = 1.0 - pow(total_angular_momentum / ((current->mass + mass) * sqrt(new_a)), 2.0);
            if (new_ecc < 0.0 || new_ecc >= 1.0)
                new_ecc = 0.0;

            double combined_mass = current->mass + mass;

            // Update current planet
            current->a = new_a;
            current->e = sqrt(new_ecc);
            current->mass = combined_mass;

            accrete_dust(&combined_mass, new_a, current->e,
                         stellar_luminosity_ratio, body_inner_bound, body_outer_bound);

            merged = TRUE;
            break;
        }

        previous = current;
        current = current->next_planet;
    }

    if (!merged)
    {
        planet_pointer new_planet = (planets *)malloc(sizeof(planets));
        new_planet->a = a;
        new_planet->e = e;
        new_planet->mass = mass;
        new_planet->gas_giant = (mass >= crit_mass);
        new_planet->next_planet = NULL;

        // Insert in sorted order
        if (planet_head == NULL || a < planet_head->a)
        {
            new_planet->next_planet = planet_head;
            planet_head = new_planet;
        }
        else
        {
            planet_pointer iterator = planet_head;
            while (iterator->next_planet != NULL && iterator->next_planet->a < a)
            {
                iterator = iterator->next_planet;
            }
            new_planet->next_planet = iterator->next_planet;
            iterator->next_planet = new_planet;
        }
    }
}



planet_pointer distribute_planetary_masses(stellar_mass_ratio,
					   stellar_luminosity_ratio, inner_dust, outer_dust)
double stellar_mass_ratio, stellar_luminosity_ratio, inner_dust, outer_dust;
{
     double a, e, mass, crit_mass,
     planetesimal_inner_bound, planetesimal_outer_bound;

     set_initial_conditions(inner_dust,outer_dust);
     planetesimal_inner_bound = innermost_planet(stellar_mass_ratio);
     planetesimal_outer_bound = outermost_planet(stellar_mass_ratio);
     while (dust_left)
     {
	  a = random_number(planetesimal_inner_bound,planetesimal_outer_bound);
	  e = random_eccentricity( );
	  mass = PROTOPLANET_MASS;
#ifdef VERBOSE
	  printf("Checking %f AU.\n",a);
#endif
	  if (dust_available(inner_effect_limit(a, e, mass),
			     outer_effect_limit(a, e, mass))) {
#ifdef VERBOSE
		    printf(".. Injecting protoplanet.\n");
#endif
		    dust_density = DUST_DENSITY_COEFF * sqrt(stellar_mass_ratio)
			 * exp(-ALPHA * pow(a,(1.0 / N)));
		    crit_mass = critical_limit(a,e,stellar_luminosity_ratio);
		    accrete_dust(&(mass),a,e,crit_mass,
				 planetesimal_inner_bound,
				 planetesimal_outer_bound);
		    if ((mass != 0.0) && (mass != PROTOPLANET_MASS))
			 coalesce_planetesimals(a,e,mass,crit_mass,
						stellar_luminosity_ratio,
						planetesimal_inner_bound,planetesimal_outer_bound);
#ifdef VERBOSE
		    else printf(".. failed due to large neighbor.\n");
#endif
	       }
#ifdef VERBOSE
	  else printf(".. failed.\n");
#endif
     }


     return(planet_head);
}

void apply_migration(planet_pointer planet)
{

    if(use_migration==1)
    {
 
        planet->original_a = planet->a;  // original orbit
        if (use_coeff_migration)
            {
                // exponential migration , 
                double decay = migration_coeff;
                planet->migrated_a = planet->original_a * decay;
            }
        if (use_exp_migration)
            {
                // exponential migration , 
                double decay = exp(migration_exp);
                planet->migrated_a = planet->original_a * decay;
            }
        if(use_random_migration)
            {
                //random migration
                double decay= (1+random_number(0,1)*0.5)*migration_coeff;
                planet->migrated_a = planet->original_a *decay;
            }
        planet->a = planet->migrated_a;
        planet->has_migrated = 1;


    }

}



planet_pointer distribute_moon_masses(planetary_mass, stellar_luminosity_ratio,
				      planet_eccentricity, inner_dust, outer_dust)
double planetary_mass, stellar_luminosity_ratio, planet_eccentricity,
     inner_dust, outer_dust;
{
     double a, e, mass, crit_mass,
     planetesimal_inner_bound, planetesimal_outer_bound;

     return(NULL);
}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// oletetaan, että nämä on määritelty jossain muualla
extern double stellar_mass_ratio, stellar_luminosity_ratio, main_seq_life, age, r_ecosphere;
extern planet_pointer planet_head;






double luminosity(mass_ratio)
double mass_ratio; 
{
     double n; 
     
     if (mass_ratio < 1.0)
	  n = 1.75 * (mass_ratio - 0.1) + 3.325;
     else 
	  n = 0.5 * (2.0 - mass_ratio) + 4.4;
     return(pow(mass_ratio,n));
}


/*--------------------------------------------------------------------------*/
/*   This function, given the orbital radius of a planet in AU, returns     */
/*   the orbital 'zone' of the particle.                                    */
/*--------------------------------------------------------------------------*/

int orbital_zone(orbital_radius)
double orbital_radius; 
{
     if (orbital_radius < (4.0 * sqrt(stellar_luminosity_ratio)))
	  return(1);
     else 
     {
	  if ((orbital_radius >= (4.0 * sqrt(stellar_luminosity_ratio))) && (orbital_radius < (15.0 * sqrt(stellar_luminosity_ratio))))
	       return(2);
	  else 
	       return(3);
     }
}


/*--------------------------------------------------------------------------*/
/*   The mass is in units of solar masses, and the density is in units      */
/*   of grams/cc.  The radius returned is in units of km.                   */
/*--------------------------------------------------------------------------*/

double volume_radius(mass, density)
double mass, density;
{
     double volume; 
     
     mass = mass * SOLAR_MASS_IN_GRAMS;
     volume = mass / density;
     return(pow((3.0 * volume) / (4.0 * PI),(1.0 / 3.0)) / CM_PER_KM);
}

/*--------------------------------------------------------------------------*/
/*    Returns the radius of the planet in kilometers.                       */
/*   The mass passed in is in units of solar masses, the orbital radius     */
/*   in A.U.                                                                */
/*   This formula is listed as eq.9 in Fogg's article, although some typos  */
/*   crop up in that eq.  See "The Internal Constitution of Planets", by    */
/*   Dr. D. S. Kothari, Mon. Not. of the Royal Astronomical Society, vol 96 */
/*   pp.833-843, 1936 for the derivation.  Specifically, this is Kothari's  */
/*   eq.23, which appears on page 840.                                      */
/*--------------------------------------------------------------------------*/

double kothari_radius(mass, orbital_radius, giant, zone)
double mass, orbital_radius;
int giant, zone;
{
     double temp, temp2, atomic_weight, atomic_num;
     
     if (zone == 1)
     {
	  if (giant)
	  {
	       atomic_weight = 9.5;
	       atomic_num = 4.5;
	  }
	  else 
	  {
	       atomic_weight = 15.0;
	       atomic_num = 8.0;
	  }
     }
     else 
	  if (zone == 2)
	  {
	       if (giant)
	       {
		    atomic_weight = 2.47;
		    atomic_num = 2.0;
	       }
	       else 
	       {
		    atomic_weight = 10.0;
		    atomic_num = 5.0;
	       }
	  }
	  else 
	  {
	       if (giant)
	       {
		    atomic_weight = 7.0;
		    atomic_num = 4.0;
	       }
	       else 
	       {
		    atomic_weight = 10.0;
		    atomic_num = 5.0;
	       }
	  }
     temp = atomic_weight * atomic_num;
     temp = (2.0 * BETA_20 * pow(SOLAR_MASS_IN_GRAMS,(1.0 / 3.0))) / (A1_20 * pow(temp,(1.0 / 3.0)));
     temp2 = A2_20 * pow(atomic_weight,(4.0 / 3.0)) * pow(SOLAR_MASS_IN_GRAMS,(2.0 / 3.0));
     temp2 = temp2 * pow(mass,(2.0 / 3.0));
     temp2 = temp2 / (A1_20 * pow(atomic_num, 2.0));
     temp2 = 1.0 + temp2;
     temp = temp / temp2;
     temp = (temp * pow(mass,(1.0 / 3.0))) / CM_PER_KM;
     return(temp);
}


/*--------------------------------------------------------------------------*/
/*  The mass passed in is in units of solar masses, and the orbital radius  */
/*  is in units of AU.  The density is returned in units of grams/cc.       */
/*--------------------------------------------------------------------------*/

double empirical_density(mass, orbital_radius, gas_giant)
double mass, orbital_radius;
int gas_giant;
{
     double temp; 
     
     temp = pow(mass * EARTH_MASSES_PER_SOLAR_MASS,(1.0 / 8.0));
     temp = temp * pow(r_ecosphere / orbital_radius,(1.0 / 4.0));
     if (gas_giant)
	  return(temp * 1.2);
     else 
	  return(temp * 5.5);
}


/*--------------------------------------------------------------------------*/
/*  The mass passed in is in units of solar masses, and the equatorial      */
/*  radius is in km.  The density is returned in units of grams/cc.         */
/*--------------------------------------------------------------------------*/

double volume_density(mass, equatorial_radius)
double mass, equatorial_radius;
{
     double volume; 
     
     mass = mass * SOLAR_MASS_IN_GRAMS;
     equatorial_radius = equatorial_radius * CM_PER_KM;
     volume = (4.0 * PI * pow(equatorial_radius,3.0)) / 3.0;
     return(mass / volume);
}


/*--------------------------------------------------------------------------*/
/*  The separation is in units of AU, and both masses are in units of solar */
/*  masses.  The period returned is in terms of Earth days.                 */
/*--------------------------------------------------------------------------*/

double period(separation, small_mass, large_mass)
double separation, small_mass, large_mass;
{
     double period_in_years; 
     
     period_in_years = sqrt(pow(separation,3.0) / (small_mass + large_mass));
     return(period_in_years * DAYS_IN_A_YEAR);
}


/*--------------------------------------------------------------------------*/
/*   Fogg's information for this routine came from Dole "Habitable Planets  */
/* for Man", Blaisdell Publishing Company, NY, 1964.  From this, he came    */
/* up with his eq.12, which is the equation for the base_angular_velocity   */
/* below.  Going a bit further, he found an equation for the change in      */
/* angular velocity per time (dw/dt) from P. Goldreich and S. Soter's paper */
/* "Q in the Solar System" in Icarus, vol 5, pp.375-389 (1966).  Comparing  */
/* to the change in angular velocity for the Earth, we can come up with an  */
/* approximation for our new planet (his eq.13) and take that into account. */
/*--------------------------------------------------------------------------*/

double day_length(mass, radius, orbital_period, eccentricity, giant)
double mass, radius, orbital_period, eccentricity;
int giant;
{
     double base_angular_velocity, planetary_mass_in_grams, k2, temp,
     equatorial_radius_in_cm, change_in_angular_velocity, spin_resonance_period;
     
     spin_resonance = FALSE;
     if (giant)
	  k2 = 0.24;
     else 
	  k2 = 0.33;
     planetary_mass_in_grams = mass * SOLAR_MASS_IN_GRAMS;
     equatorial_radius_in_cm = radius * CM_PER_KM;
     base_angular_velocity = sqrt(2.0 * J * (planetary_mass_in_grams) / (k2 * pow(equatorial_radius_in_cm, 2.0)));
     /*   This next term describes how much a planet's rotation is slowed by    */
     /*  it's moons.  Find out what dw/dt is after figuring out Goldreich and   */
     /*  Soter's Q'.                                                            */
     change_in_angular_velocity = 0.0;
     temp = base_angular_velocity + (change_in_angular_velocity * age);
     /*   'temp' is now the angular velocity. Now we change from rad/sec to     */
     /*  hours/rotation.							   */
     temp = 1.0 / ((temp / radians_per_rotation) * SECONDS_PER_HOUR);
     if (temp >= orbital_period)
     {
	  spin_resonance_period = ((1.0 - eccentricity) / (1.0 + eccentricity)) * orbital_period;
	  if (eccentricity > 0.1)
	  {
	       temp = spin_resonance_period;
	       spin_resonance = TRUE;
	  }
	  else 
	       temp = orbital_period;
     }
     return(temp);
}


/*--------------------------------------------------------------------------*/
/*   The orbital radius is expected in units of Astronomical Units (AU).    */
/*   Inclination is returned in units of degrees.                           */
/*--------------------------------------------------------------------------*/

int inclination(orbital_radius)
double orbital_radius; 
{
     int temp; 
     
     temp = (int)(pow(orbital_radius,0.2) * about(EARTH_AXIAL_TILT,0.4));
     return(temp % 360);
}


/*--------------------------------------------------------------------------*/
/*   This function implements the escape velocity calculation.  Note that   */
/*  it appears that Fogg's eq.15 is incorrect.                              */
/*  The mass is in units of solar mass, the radius in kilometers, and the   */
/*  velocity returned is in cm/sec.                                         */
/*--------------------------------------------------------------------------*/

double escape_vel(mass, radius)
double mass, radius;
{
     double mass_in_grams, radius_in_cm;
     
     mass_in_grams = mass * SOLAR_MASS_IN_GRAMS;
     radius_in_cm = radius * CM_PER_KM;
     return(sqrt(2.0 * GRAV_CONSTANT * mass_in_grams / radius_in_cm));
}


/*--------------------------------------------------------------------------*/
/*  This is Fogg's eq.16.  The molecular weight (usually assumed to be N2)  */
/*  is used as the basis of the Root Mean Square velocity of the molecule   */
/*  or atom.  The velocity returned is in cm/sec.                           */
/*--------------------------------------------------------------------------*/

double rms_vel(molecular_weight, orbital_radius)
double molecular_weight, orbital_radius;
{
     double exospheric_temp; 
     
     exospheric_temp = EARTH_EXOSPHERE_TEMP / pow(orbital_radius, 2.0);
     return(sqrt((3.0 * MOLAR_GAS_CONST * exospheric_temp) / molecular_weight) * CM_PER_METER);
}


/*--------------------------------------------------------------------------*/
/*   This function returns the smallest molecular weight retained by the    */
/*  body, which is useful for determining the atmosphere composition.       */
/*  Orbital radius is in A.U.(ie: in units of the earth's orbital radius),  *)
    (*  mass is in units of solar masses, and equatorial radius is in units of  */
/*  kilometers.                                                             */
/*--------------------------------------------------------------------------*/

double molecule_limit(orbital_radius, mass, equatorial_radius)
double orbital_radius, mass, equatorial_radius;
{
     double numerator, denominator1, denominator2, escape_velocity, temp;
     
     escape_velocity = escape_vel(mass,equatorial_radius);
     return((3.0 * pow(GAS_RETENTION_THRESHOLD * CM_PER_METER, 2.0) * MOLAR_GAS_CONST * EARTH_EXOSPHERE_TEMP) / pow(escape_velocity, 2.0));
}


/*--------------------------------------------------------------------------*/
/*   This function calculates the surface acceleration of a planet.  The    */
/*  mass is in units of solar masses, the radius in terms of km, and the    */
/*  acceleration is returned in units of cm/sec2.                           */
/*--------------------------------------------------------------------------*/

double acceleration(mass, radius)
double mass, radius;
{
     return(GRAV_CONSTANT * (mass * SOLAR_MASS_IN_GRAMS) / pow(radius * CM_PER_KM, 2.0));
}


/*--------------------------------------------------------------------------*/
/*   This function calculates the surface gravity of a planet.  The         */
/*  acceleration is in units of cm/sec2, and the gravity is returned in     */
/*  units of Earth gravities.                                               */
/*--------------------------------------------------------------------------*/

double gravity(acceleration)
double acceleration; 
{
     return(acceleration / EARTH_ACCELERATION);
}


/*--------------------------------------------------------------------------*/
/*  Note that if the orbital radius of the planet is greater than or equal  */
/*  to R_inner, 99% of it's volatiles are assumed to have been deposited in */
/*  surface reservoirs (otherwise, it suffers from the greenhouse effect).  */
/*--------------------------------------------------------------------------*/

int greenhouse(zone, orbital_radius, greenhouse_radius)
int zone; 
double orbital_radius, greenhouse_radius;
{
     if ((orbital_radius < greenhouse_radius) && (zone == 1))
	  return(TRUE);
     else 
	  return(FALSE);
}


/*--------------------------------------------------------------------------*/
/*  This implements Fogg's eq.17.  The 'inventory' returned is unitless.    */
/*--------------------------------------------------------------------------*/

double vol_inventory(mass, escape_vel, rms_vel, stellar_mass, zone, greenhouse_effect)
double mass, escape_vel, rms_vel, stellar_mass;
int zone, greenhouse_effect;
{
     double velocity_ratio, proportion_const, temp1, temp2, mass_in_earth_units;
     
     velocity_ratio = escape_vel / rms_vel;
     if (velocity_ratio >= GAS_RETENTION_THRESHOLD)
     {
	  switch (zone) {
	       case 1:
		    proportion_const = 100000.0;
		    break;
	       case 2:
		    proportion_const = 75000.0;
		    break;
	       case 3:
		    proportion_const = 250.0;
		    break;
	       default:
		    printf("Error: orbital zone not initialized correctly!\n");
		    break;
	       }
	  mass_in_earth_units = mass * EARTH_MASSES_PER_SOLAR_MASS;
	  temp1 = (proportion_const * mass_in_earth_units) / stellar_mass;
	  temp2 = about(temp1,0.2);
	  if (greenhouse_effect)
	       return(temp2);
	  else 
	       return(temp2 / 100.0);
     }
     else 
	  return(0.0);
}


/*--------------------------------------------------------------------------*/
/*  This implements Fogg's eq.18.  The pressure returned is in units of     */
/*  millibars (mb).  The gravity is in units of Earth gravities, the radius */
/*  in units of kilometers.                                                 */
/*--------------------------------------------------------------------------*/

double pressure(volatile_gas_inventory, equatorial_radius, gravity)
double volatile_gas_inventory, equatorial_radius, gravity;
{
     equatorial_radius = EARTH_RADIUS_IN_KM / equatorial_radius;
     return(volatile_gas_inventory * gravity / pow(equatorial_radius, 2.0));
}

/*--------------------------------------------------------------------------*/
/*   This function returns the boiling point of water in an atmosphere of   */
/*   pressure 'surface_pressure', given in millibars.  The boiling point is */
/*   returned in units of Kelvin.  This is Fogg's eq.21.                    */
/*--------------------------------------------------------------------------*/

double boiling_point(surface_pressure)
double surface_pressure; 
{
     double surface_pressure_in_bars; 
     
     surface_pressure_in_bars = surface_pressure / MILLIBARS_PER_BAR;
     return(1.0 / (log(surface_pressure_in_bars) / -5050.5 + 1.0 / 373.0));
}


/*--------------------------------------------------------------------------*/
/*   This function is Fogg's eq.22.  Given the volatile gas inventory and   */
/*   planetary radius of a planet (in Km), this function returns the        */
/*   fraction of the planet covered with water.                             */
/*   I have changed the function very slightly:  the fraction of Earth's    */
/*   surface covered by water is 71%, not 75% as Fogg used.                 */
/*--------------------------------------------------------------------------*/

double hydrosphere_fraction(volatile_gas_inventory, planetary_radius)
double volatile_gas_inventory, planetary_radius;
{
     double temp; 
     
     temp = (0.71 * volatile_gas_inventory / 1000.0) * pow(EARTH_RADIUS_IN_KM / planetary_radius, 2.0);
     if (temp >= 1.0)
	  return(1.0);
     else 
	  return(temp);
}


/*--------------------------------------------------------------------------*/
/*   Given the surface temperature of a planet (in Kelvin), this function   */
/*   returns the fraction of cloud cover available.  This is Fogg's eq.23.  */
/*   See Hart in "Icarus" (vol 33, pp23 - 39, 1978) for an explanation.     */
/*   This equation is Hart's eq.3.                                          */
/*   I have modified it slightly using constants and relationships from     */
/*   Glass's book "Introduction to Planetary Geology", p.46.                */
/*   The 'CLOUD_COVERAGE_FACTOR' is the amount of surface area on Earth     */
/*   covered by one Kg. of cloud.					    */
/*--------------------------------------------------------------------------*/

double cloud_fraction(surface_temp, smallest_MW_retained, equatorial_radius, hydrosphere_fraction)
double surface_temp, smallest_MW_retained, equatorial_radius,
     hydrosphere_fraction;
{
     double water_vapor_in_kg, fraction, surface_area, hydrosphere_mass;
     
     if (smallest_MW_retained > WATER_VAPOR)
	  return(0.0);
     else 
     {
	  surface_area = 4.0 * PI * pow(equatorial_radius, 2.0);
	  hydrosphere_mass = hydrosphere_fraction * surface_area * EARTH_WATER_MASS_PER_AREA;
	  water_vapor_in_kg = (0.00000001 * hydrosphere_mass) * exp(Q2_36 * (surface_temp - 288.0));
	  fraction = CLOUD_COVERAGE_FACTOR * water_vapor_in_kg / surface_area;
	  if (fraction >= 1.0)
	       return(1.0);
	  else 
	       return(fraction);
     }
}


/*--------------------------------------------------------------------------*/
/*   Given the surface temperature of a planet (in Kelvin), this function   */
/*   returns the fraction of the planet's surface covered by ice.  This is  */
/*   Fogg's eq.24.  See Hart[24] in Icarus vol.33, p.28 for an explanation. */
/*   I have changed a constant from 70 to 90 in order to bring it more in   */
/*   line with the fraction of the Earth's surface covered with ice, which  */
/*   is approximatly .016 (=1.6%).                                          */
/*--------------------------------------------------------------------------*/

double ice_fraction(hydrosphere_fraction, surface_temp)
double hydrosphere_fraction, surface_temp;
{
     double temp; 
     
     if (surface_temp > 328.0) 
	  surface_temp = 328.0;
     temp = pow(((328.0 - surface_temp) / 90.0),5.0);
     if (temp > (1.5 * hydrosphere_fraction))
	  temp = (1.5 * hydrosphere_fraction);
     if (temp >= 1.0)
	  return(1.0);
     else 
	  return(temp);
}


/*--------------------------------------------------------------------------*/
/*  This is Fogg's eq.19.  The ecosphere radius is given in AU, the orbital */
/*  radius in AU, and the temperature returned is in Kelvin.		    */
/*--------------------------------------------------------------------------*/

double eff_temp(ecosphere_radius, orbital_radius, albedo)
double ecosphere_radius, orbital_radius, albedo;
{
     return(sqrt(ecosphere_radius / orbital_radius) * pow((1.0 - albedo) / 0.7,0.25) * EARTH_EFFECTIVE_TEMP);
}


/*--------------------------------------------------------------------------*/
/*  This is Fogg's eq.20, and is also Hart's eq.20 in his "Evolution of     */
/*  Earth's Atmosphere" article.  The effective temperature given is in     */
/*  units of Kelvin, as is the rise in temperature produced by the          */
/*  greenhouse effect, which is returned.                                   */
/*--------------------------------------------------------------------------*/

double green_rise(optical_depth, effective_temp, surface_pressure)
double optical_depth, effective_temp, surface_pressure;
{
     double convection_factor; 
     
     convection_factor = EARTH_CONVECTION_FACTOR * pow((surface_pressure / EARTH_SURF_PRES_IN_MILLIBARS),0.25);
     return(pow((1.0 + 0.75 * optical_depth),0.25) - 1.0) * effective_temp * convection_factor;
}


/*--------------------------------------------------------------------------*/
/*   The surface temperature passed in is in units of Kelvin.               */
/*   The cloud adjustment is the fraction of cloud cover obscuring each     */
/*   of the three major components of albedo that lie below the clouds.     */
/*--------------------------------------------------------------------------*/

double planet_albedo(water_fraction, cloud_fraction, ice_fraction, surface_pressure)
double water_fraction, cloud_fraction, ice_fraction, surface_pressure;
{
     double rock_fraction, cloud_adjustment, components, cloud_contribution,
     rock_contribution, water_contribution, ice_contribution;
     
     rock_fraction = 1.0 - water_fraction - ice_fraction;
     components = 0.0;
     if (water_fraction > 0.0)
	  components = components + 1.0;
     if (ice_fraction > 0.0)
	  components = components + 1.0;
     if (rock_fraction > 0.0)
	  components = components + 1.0;
     cloud_adjustment = cloud_fraction / components;
     if (rock_fraction >= cloud_adjustment)
	  rock_fraction = rock_fraction - cloud_adjustment;
     else 
	  rock_fraction = 0.0;
     if (water_fraction > cloud_adjustment)
	  water_fraction = water_fraction - cloud_adjustment;
     else 
	  water_fraction = 0.0;
     if (ice_fraction > cloud_adjustment)
	  ice_fraction = ice_fraction - cloud_adjustment;
     else 
	  ice_fraction = 0.0;
     cloud_contribution = cloud_fraction * about(CLOUD_ALBEDO,0.2);
     if (surface_pressure == 0.0)
	  rock_contribution = rock_fraction * about(AIRLESS_ROCKY_ALBEDO,0.3);
     else 
	  rock_contribution = rock_fraction * about(ROCKY_ALBEDO,0.1);
     water_contribution = water_fraction * about(WATER_ALBEDO,0.2);
     if (surface_pressure == 0.0)
	  ice_contribution = ice_fraction * about(AIRLESS_ICE_ALBEDO,0.4);
     else 
	  ice_contribution = ice_fraction * about(ICE_ALBEDO,0.1);
     return(cloud_contribution + rock_contribution + water_contribution + ice_contribution);
}


/*--------------------------------------------------------------------------*/
/*   This function returns the dimensionless quantity of optical depth,     */
/*   which is useful in determining the amount of greenhouse effect on a    */
/*   planet.                                                                */
/*--------------------------------------------------------------------------*/

double opacity(molecular_weight, surface_pressure)
double molecular_weight, surface_pressure;
{
     double optical_depth; 
     
     optical_depth = 0.0;
     if ((molecular_weight >= 0.0) && (molecular_weight < 10.0))
	  optical_depth = optical_depth + 3.0;
     if ((molecular_weight >= 10.0) && (molecular_weight < 20.0))
	  optical_depth = optical_depth + 2.34;
     if ((molecular_weight >= 20.0) && (molecular_weight < 30.0))
	  optical_depth = optical_depth + 1.0;
     if ((molecular_weight >= 30.0) && (molecular_weight < 45.0))
	  optical_depth = optical_depth + 0.15;
     if ((molecular_weight >= 45.0) && (molecular_weight < 100.0))
	  optical_depth = optical_depth + 0.05;
     if (surface_pressure >= (70.0 * EARTH_SURF_PRES_IN_MILLIBARS))
	  optical_depth = optical_depth * 8.333;
     else 
	  if (surface_pressure >= (50.0 * EARTH_SURF_PRES_IN_MILLIBARS))
	       optical_depth = optical_depth * 6.666;
	  else 
	       if (surface_pressure >= (30.0 * EARTH_SURF_PRES_IN_MILLIBARS))
		    optical_depth = optical_depth * 3.333;
	       else 
		    if (surface_pressure >= (10.0 * EARTH_SURF_PRES_IN_MILLIBARS))
			 optical_depth = optical_depth * 2.0;
		    else 
			 if (surface_pressure >= (5.0 * EARTH_SURF_PRES_IN_MILLIBARS))
			      optical_depth = optical_depth * 1.5;
     return(optical_depth);
}


/*--------------------------------------------------------------------------*/
/*   The temperature calculated is in degrees Kelvin.                       */
/*   Quantities already known which are used in these calculations:         */
/*	 planet->molecule_weight					    */
/*	 planet->surface_pressure					    */
/*       R_ecosphere                                                        */
/*	 planet->a							    */
/*	 planet->volatile_gas_inventory					    */
/*	 planet->radius							    */
/*	 planet->boil_point						    */
/*--------------------------------------------------------------------------*/

void iterate_surface_temp(planet)
planet_pointer *planet; 
{
     double surface_temp, effective_temp, greenhouse_rise, previous_temp,
     optical_depth, albedo, water, clouds, ice;
     
     optical_depth = opacity((*planet)->molecule_weight,(*planet)->surface_pressure);
     effective_temp = eff_temp(r_ecosphere,(*planet)->a,EARTH_ALBEDO);
     greenhouse_rise = green_rise(optical_depth,effective_temp,(*planet)->surface_pressure);
     surface_temp = effective_temp + greenhouse_rise;
     previous_temp = surface_temp - 5.0;		/* force the while loop the first time */
    int nn=0;
     while ((fabs(surface_temp - previous_temp) > 1.0)) {
	       previous_temp = surface_temp;
	       water = hydrosphere_fraction((*planet)->volatile_gas_inventory,(*planet)->radius);
	       clouds = cloud_fraction(surface_temp,(*planet)->molecule_weight,(*planet)->radius,water);
	       ice = ice_fraction(water,surface_temp);
	       if ((surface_temp >= (*planet)->boil_point) || (surface_temp <= FREEZING_POINT_OF_WATER))
		    water = 0.0;
	       albedo = planet_albedo(water,clouds,ice,(*planet)->surface_pressure);
	       optical_depth = opacity((*planet)->molecule_weight,(*planet)->surface_pressure);
	       effective_temp = eff_temp(r_ecosphere,(*planet)->a,albedo);
	       greenhouse_rise = green_rise(optical_depth,effective_temp,(*planet)->surface_pressure);
	       surface_temp = effective_temp + greenhouse_rise;
            nn++;
           if(nn>200) 
            {
            // effective_temp = eff_temp(r_ecosphere,(*planet)->a,albedo);
            // surface_temp = effective_temp ;
                break;
            }
	  }
     (*planet)->hydrosphere = water;
     (*planet)->cloud_cover = clouds;
     (*planet)->ice_cover = ice;
     (*planet)->albedo = albedo;
     (*planet)->surface_temp = surface_temp;
}
#include	<stdio.h>
#include	<string.h>
#include        <math.h>

#ifdef MSDOS
#include        <stddef.h>
#include        <malloc.h>
#include	<stdlib.h>
#include        <float.h>
#endif



/*#define	VERBOSE*/

/*  These are all of the global variables used during accretion:  */
planet_pointer planet_head;
double stellar_mass_ratio, stellar_luminosity_ratio, main_seq_life,
     age, r_ecosphere, r_greenhouse, radians_per_rotation;
int spin_resonance;




void init()
{
  //   srand(25);
}

void generate_stellar_system()
{
     planet_pointer planet; 

     radians_per_rotation = 2.0 * PI;
     stellar_mass_ratio = random_number(0.6,1.3);
     stellar_luminosity_ratio = luminosity(stellar_mass_ratio);
     planet = distribute_planetary_masses(stellar_mass_ratio,stellar_luminosity_ratio,0.0,stellar_dust_limit(stellar_mass_ratio));
     main_seq_life = 1.0E10 * (stellar_mass_ratio / stellar_luminosity_ratio);
     if ((main_seq_life >= 6.0E9))
	  age = random_number(1.0E9,6.0E9);
     else 
	  age = random_number(1.0E9,main_seq_life);
     r_ecosphere = sqrt(stellar_luminosity_ratio);
     r_greenhouse = r_ecosphere * GREENHOUSE_EFFECT_CONST;
     while (planet != NULL)
     {
     
/*	planet->first_moon = distribute_moon_masses (planet->mass,
		stellar_luminosity_ratio, planet->e,
		0.0, planet_dust_limit(planet->mass));*/
	  planet->orbit_zone = orbital_zone(planet->a);
	  if (planet->gas_giant)
	  {
        
	       planet->density = empirical_density(planet->mass,planet->a,planet->gas_giant);
	       planet->radius = volume_radius(planet->mass,planet->density);
	  }
	  else 
	  {
	       planet->radius = kothari_radius(planet->mass,planet->a,planet->gas_giant,planet->orbit_zone);
	       planet->density = volume_density(planet->mass,planet->radius);
	  }
	  
 // apply_migration(planet);
planet->orbital_period = period(planet->a,planet->mass,stellar_mass_ratio);
	  planet->day = day_length(planet->mass,planet->radius,planet->orbital_period,planet->e,planet->gas_giant);
	  planet->resonant_period = spin_resonance;
	  planet->axial_tilt = inclination(planet->a);
	  planet->escape_velocity = escape_vel(planet->mass,planet->radius);
	  planet->surface_accel = acceleration(planet->mass,planet->radius);
	  planet->rms_velocity = rms_vel(MOLECULAR_NITROGEN,planet->a);
	  planet->molecule_weight = molecule_limit(planet->a,planet->mass,planet->radius);
	  if ((planet->gas_giant))
	  {
	       planet->surface_grav = INCREDIBLY_LARGE_NUMBER;
	       planet->greenhouse_effect = FALSE;
	       planet->volatile_gas_inventory = INCREDIBLY_LARGE_NUMBER;
	       planet->surface_pressure = INCREDIBLY_LARGE_NUMBER;
	       planet->boil_point = INCREDIBLY_LARGE_NUMBER;
	       planet->hydrosphere = INCREDIBLY_LARGE_NUMBER;
	       planet->albedo = about(GAS_GIANT_ALBEDO,0.1);
	       planet->surface_temp = INCREDIBLY_LARGE_NUMBER;
	  }
	  else 
	  {
	       planet->surface_grav = gravity(planet->surface_accel);
	       planet->greenhouse_effect = greenhouse(planet->orbit_zone,planet->a,r_greenhouse);
	       planet->volatile_gas_inventory = vol_inventory(planet->mass,planet->escape_velocity,planet->rms_velocity,stellar_mass_ratio,planet->orbit_zone,planet->greenhouse_effect);
	       planet->surface_pressure = pressure(planet->volatile_gas_inventory,planet->radius,planet->surface_grav);
	       if ((planet->surface_pressure == 0.0))
		    planet->boil_point = 0.0;
	       else 
		    planet->boil_point = boiling_point(planet->surface_pressure);
	       iterate_surface_temp(&(planet));
	  }

     apply_migration(planet);
    
  iterate_surface_temp(&(planet));

	  planet = planet->next_planet;

   
//{  
 
    }
     display_system( );
}

/*
int main (int argc,char *argv[])

{
     init( );
     generate_stellar_system( );
return(1);

}
*/



int main(int argc, char *argv[]) {
    // Set default values for your simulation parameters
  
    printf("Program name: %s\n", argv[0]);

   

    // Loop through command-line arguments, starting from the second one (argv[0] is program name)
    for (int i = 1; i < argc; i++) {
        // Check for help options (-h or --help)
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("\nUsage:\n");
            printf("  %s [-mass <value>] [-migrate] [-coeff <value>] [-ecosphere <value>]\n", argv[0]);
            printf("  -mass <value>      : Stellar mass ratio (double, default 1.0)\n");
            printf("  -migrate           : Enable planetary migration (no argument needed)\n");
            printf("  -coeff <value>     : Migration coefficient (double, default 0.1)\n");
            printf("  -ecosphere <value> : Ecosphere radius (double, default 1.0)\n");
            printf("  -h, --help         : Display this help message\n");
            return 0; // Exit after showing help
        }
        // Handle -mass parameter
        else if (strcmp(argv[i], "-mass") == 0) {
            if (i + 1 < argc) { // Ensure there's a value for the argument
                stellar_mass_ratio = atof(argv[i+1]); // Convert string to double
                printf("Stellar mass ratio set to: %.2f\n", stellar_mass_ratio);
                i++; // Increment i to skip the value we just processed (argv[i+1])
            } else {
                fprintf(stderr, "Error: -mass requires a value.\n");
                return 1; // Indicate an error
            }
        }
        // Handle -migrate (flag) parameter
        else if (strcmp(argv[i], "-migrate") == 0) {
            use_migration = 1;
            printf("Planetary migration enabled.\n");
        }
        // Handle -coeff parameter
        else if (strcmp(argv[i], "-coeff") == 0) {
            if (i + 1 < argc) {
                migration_coeff = atof(argv[i+1]);
                printf("Migration coefficient set to: %.2f\n", migration_coeff);
                i++;
            } else {
                fprintf(stderr, "Error: -coeff requires a value.\n");
                return 1;
            }
        }
        // Handle -ecosphere parameter
        else if (strcmp(argv[i], "-ecosphere") == 0) {
            if (i + 1 < argc) {
                r_ecosphere = atof(argv[i+1]);
                printf("Ecosphere radius set to: %.2f\n", r_ecosphere);
                i++;
            } else {
                fprintf(stderr, "Error: -ecosphere requires a value.\n");
                return 1;
            }
        }
        else if (strcmp(argv[i], "-seed") == 0) {
            if (i + 1 < argc) {
                seed1= atoi(argv[i+1]);
                     srand(seed1);
                printf(" Seed set to: %.2f\n", r_ecosphere);
                i++;
            } else {
                fprintf(stderr, "Random seed\n");
                srand(time(NULL)); 
                //randomize();
                return 1;
            }
        } 
/*
double base_dust_limit=200;
double base_innermost_planet=0.3;
*/  
        else if (strcmp(argv[i], "-disk_radius") == 0) {
            if (i + 1 < argc) {
                base_dust_limit= atof(argv[i+1]);
         
                printf("Disk radius set to: %.2f\n",  base_dust_limit);
                i++;
            } else {
                fprintf(stderr, " -disk_radius requires a parameter value.\n");
       
                return 1;
            }
        } 

        else if (strcmp(argv[i], "-innermost_planet") == 0) {
            if (i + 1 < argc) {
                 base_innermost_planet= atof(argv[i+1]);
         
                printf("innermost planet set to: %.2f\n",   base_innermost_planet);
                i++;
            } else {
                fprintf(stderr, "-innermost_planet equires a parameter value.\n");
       
                return 1;
            }
        } 

     else {
            fprintf(stderr, "Unknown parameter: %s. Use -h for help.\n", argv[i]);
            return 1; // Indicate an error
        }
    }

    printf("\nUsed parameter values:\n");
    printf("  Random number seed:    %i   \n", seed1);
    printf("  Stellar Mass Ratio:    %.2f \n", stellar_mass_ratio);
    printf("  Innermost planet:      %i   \n",  base_innermost_planet );
    printf("  Disk radius :          %i   \n",  base_dust_limit );

    printf("  Use Migration:         %d\n", use_migration);
    printf("  Migration coefficient: %.2f\n", migration_coeff);
    printf("  Ecosphere radius:      %.2f\n", r_ecosphere);
   

      init( );
     generate_stellar_system( );

        if(render==1)
        {
            generate_and_render_povray() ;
      }   
 return 0; // Program executed successfully
}

int display_system(void)
{
     planet_pointer node1;
     int counter;
     FILE *fp;

     printf("                         SYSTEM  CHARACTERISTICS\n");
     printf("Mass of central star (in solar masses): %4.2lf\n", stellar_mass_ratio);
     printf("Luminosity of central star (relative to the sun): %5.2lf\n", stellar_luminosity_ratio);
     printf("Total main sequence lifetime (in million yrs): %10.3lf\n", (main_seq_life / 1.0E6));
     printf("Current age of stellar system (in million yrs): %10.3lf\n", (age / 1.0E6));
     printf("Radius of habitable ecosphere (AU): %3.3lf\n", r_ecosphere);

     // Avaa tiedosto ja kirjoita otsikot
     fp = fopen("system_output.csv", "w");
     if (fp == NULL) {
         perror("Unable to open CSV file");
         return -1;
     }

     fprintf(fp, "Planet,Distance_AU,Eccentricity,Mass_Earth,Radius_km,Density_gcc,Escape_kms,"
                 "Mol_weight,Accel_cm_s2,Grav_g,Boil_C,Pressure_atm,Temp_C,Hydro_%%,Cloud_%%,"
                 "Ice_%%,Tilt_deg,Albedo,Year_days,Day_hours,Gas_giant,Resonant\n");

     node1 = planet_head;
     counter = 1;
     while (node1 != NULL)
     {
         printf("Planet #%d:\n", counter);
         if (node1->gas_giant)
             printf("Gas giant...\n");
         if (node1->resonant_period)
             printf("In resonant period with primary.\n");

         // Tulostus terminaaliin kuten aiemmin (säilytetty)
         printf("   Distance from primary star (in A.U.): %7.3lf\n", node1->a);
         printf("   Has migrated: %i \n", node1->has_migrated);  
         printf("   Eccentricity of orbit: %5.3lf\n", node1->e);
         printf("   Mass (in Earth masses): %7.3lf\n", node1->mass * EARTH_MASSES_PER_SOLAR_MASS);
         printf("   Equatorial radius (in Km): %10.1lf\n", node1->radius);
         printf("   Density (in g/cc): %6.3lf\n", node1->density);
         printf("   Escape Velocity (in km/sec): %5.2lf\n", node1->escape_velocity / CM_PER_KM);
         printf("   Smallest molecular weight retained: %5.2lf\n", node1->molecule_weight);
         printf("   Surface acceleration (in cm/sec2): %6.2lf\n", node1->surface_accel);
    if (node1->has_migrated)
    {
        printf("   *** Migrated inward from original orbit: %.3lf AU ***\n", node1->original_a);
    }
         // CSV:llä varmistetaan että kentät ovat olemassa
         double surface_grav = node1->surface_grav;
         double boil_c = node1->boil_point - KELVIN_CELCIUS_DIFFERENCE;
         double pressure_atm = node1->surface_pressure / 1000.0;
         double temp_c = node1->surface_temp - KELVIN_CELCIUS_DIFFERENCE;
         double hydro = node1->hydrosphere * 100.0;
         double cloud = node1->cloud_cover * 100.0;
         double ice = node1->ice_cover * 100.0;

         // Kirjoitus tiedostoon
         fprintf(fp, "%d,%.3lf,%.3lf,%.3lf,%.1lf,%.3lf,%.2lf,%.2lf,%.2lf,",
             counter, node1->a, node1->e,
             node1->mass * EARTH_MASSES_PER_SOLAR_MASS,
             node1->radius, node1->density,
             node1->escape_velocity / CM_PER_KM,
             node1->molecule_weight,
             node1->surface_accel);

         if (!(node1->gas_giant)) {
             fprintf(fp, "%.2lf,%.1lf,%.3lf,%.2lf,%.2lf,%.2lf,%.2lf,",
                 surface_grav, boil_c, pressure_atm,
                 temp_c, hydro, cloud, ice);
         } else {
             fprintf(fp, ",,,,,,,"); // Täytetään tyhjät kentät kaasujättiläisille
         }

         fprintf(fp, "%d,%.3lf,%.2lf,%2lf,%s,%s\n",
             node1->axial_tilt,
             node1->albedo,
             node1->orbital_period,
             node1->day,
             node1->gas_giant ? "yes" : "no",
             node1->resonant_period ? "yes" : "no");

         counter++;
         node1 = node1->next_planet;
     }

     fclose(fp);
}

/*----------------------------------------------------------------------*/
/*  This function returns a random real number between the specified    */
/* inner and outer bounds.                                              */
/*----------------------------------------------------------------------*/




double random_number(inner, outer)
double inner, outer;
{
     double range;

     range = outer - inner;
     return((((double)rand()) / (double)(RAND_MAX)) * range + inner);
}


/*----------------------------------------------------------------------*/
/*   This function returns a value within a certain variation of the    */
/*   exact value given it in 'value'.                                   */
/*----------------------------------------------------------------------*/

double about(value, variation)
double value, variation;
{
     return(value + (value * random_number(-variation,variation)));
}

double random_eccentricity()
{
     return(1.0 - pow(random_number(0.0, 1.0),ECCENTRICITY_COEFF));
}


/// muokkaus



int generate_and_render_povray() {
    FILE *fp_pov=NULL;
    char pov_filename[] = "system_scene.pov";
    char output_image_filename[] = "system_render.png";
    int counter=0;
    planet_pointer node1;

    // --- POV-Ray File Generation ---
    fp_pov = fopen(pov_filename, "w");
    if (fp_pov == NULL) {
        perror("Error: Unable to open POV-Ray file");
        return -1;
    }

    // POV-Ray header and basic scene setup
    fprintf(fp_pov, "// POV-Ray Scene file generated by System Simulator\n");
    fprintf(fp_pov, "// Generated on: %s\n", __DATE__);
    fprintf(fp_pov, "\n");
    fprintf(fp_pov, "#include \"colors.inc\"\n"); // Standard POV-Ray color definitions
    fprintf(fp_pov, "#include \"functions.inc\"\n");
    // Camera setup (adjust position and look_at for best view)
    // Positioned to look down slightly on the orbital plane
    fprintf(fp_pov, "camera {\n");
    fprintf(fp_pov, "  location <14, 0, -100> // x, y (up), z (depth)\n");
    fprintf(fp_pov, "  look_at <14, 0, 0>\n");
    fprintf(fp_pov, "  right x * image_width / image_height\n");
    fprintf(fp_pov, "  angle 20 // Field of view\n");
    fprintf(fp_pov, "}\n\n");

    // Light source - a strong white light to illuminate the scene
    // This could also be your central star itself
    fprintf(fp_pov, "light_source { <0, 0, -10> color White * 1.5 }\n\n");

    // Ground plane (optional, but makes the scene look better)
   // fprintf(fp_pov, "plane {\n");
   // fprintf(fp_pov, "  y, -0.1 // Normal vector (y-axis), distance from origin\n");
   // fprintf(fp_pov, "  pigment { checker color DarkGreen, color Green }\n");
   // fprintf(fp_pov, "}\n\n");

    // --- Central Star ---
    // Scale its radius based on stellar_mass_ratio or stellar_luminosity_ratio
    // For simplicity, let's use a fixed size that looks good with planets scaled below
    double star_radius_pov = 1.5; // POV-Ray units
    fprintf(fp_pov, "sphere {\n");
    fprintf(fp_pov, "  <0, 0, 0>, %f // Center and radius\n", star_radius_pov);
    fprintf(fp_pov, "  pigment { color <1,1,0.5> filter 0.5 } // Yellow, semi-transparent for glow effect\n");
    fprintf(fp_pov, "  // Add an emissive finish for a glowing effect\n");
    fprintf(fp_pov, "  finish { ambient 1 diffuse 0 emission 1 }\n");
    fprintf(fp_pov, "  // Add a halo for a more realistic star glow (requires photons in render settings)\n");
    // fprintf(fp_pov, "  photons { emission 1 }\n"); // Uncomment if you want to use photons for glow
    fprintf(fp_pov, "}\n\n");


    // --- Planets ---
    // Scale factor for planet radii and distances in POV-Ray units
    // Adjust these to make planets visible and to fit within the view
    double planet_radius_scale = 0.004; // 1 Km in simulation = 0.0005 POV-Ray units
    double distance_scale = 2.0;        // 1 AU in simulation = 5.0 POV-Ray units


    node1 = planet_head;
    counter = 1;
    while (node1 != NULL) {
        // Calculate planet position in POV-Ray coordinates
        // Simple circular orbit for visualization. For eccentricity, you'd need more complex math.
        // We'll place them along the X-axis for simplicity in this example.
        //double planet_x =2+ log(node1->a) * distance_scale;
      double planet_x =2+counter*distance_scale;  
      double planet_y = 0.0; // Assume they are on the x-z plane (y=0)
        double planet_z = 0.0;

        double planet_pov_radius = sqrt(node1->radius) * planet_radius_scale;
      char buffu [256]="";
   memset(buffu,0,256);

        fprintf(fp_pov, "// Planet #%d\n", counter);
        fprintf(fp_pov, "sphere {\n");
        fprintf(fp_pov, "  <%f, %f, %f>, %f // Position and radius\n",
                planet_x, planet_y, planet_z, planet_pov_radius);

        if (node1->gas_giant) {
            fprintf(fp_pov, "  pigment { color rgb <0.4, 0.8, 0.6> } // Geenish for gas giants\n");
            fprintf(fp_pov, "  finish { phong 0.8 } // Shiny finish\n");
        } else {
            // Terrestrial planet coloring based on features (simplified)
            if (node1->hydrosphere > 0.5 && node1->cloud_cover < 0.5) {
                fprintf(fp_pov, "  pigment { color rgb <0.2, 0.4, 0.8> } // Water world (blue)\n");
            } else if (node1->ice_cover > 0.3) {
                fprintf(fp_pov, "  pigment { color White } // Ice world\n");
            } else if (node1->surface_pressure > 0.001) { // Basic atmospheric planet
                fprintf(fp_pov, "  pigment { color rgb <0.2, 0.6, 0.2> } // Green/Earth-like\n");
            } else {
                fprintf(fp_pov, "  pigment { color rgb <0.7, 0.4, 0.2> } // Rocky/desert planet\n");
            }
            fprintf(fp_pov, "  finish { phong 0.6 roughness 0.05 ambient 0}\n"); // Less shiny, more diffuse
        }
        fprintf(fp_pov, "}\n\n");

       fprintf(fp_pov, "text { \n\n");
      fprintf(fp_pov, "ttf \"timrom.ttf\" ");

    memset(buffu,0,256);
  sprintf( buffu ,"\"%d\" ", counter );
  fprintf(fp_pov, buffu);
  fprintf(fp_pov, " 0.15,0 \n pigment {color rgb <1,1,1> }\n");
      fprintf(fp_pov, " translate y*-2.25 \n\n");
    fprintf(fp_pov, " translate x*%f \n\n", planet_x-0.25); // planet_x
      fprintf(fp_pov, " }\n\n");

        counter++;
        node1 = node1->next_planet;
    }

    fclose(fp_pov);
    printf("POV-Ray scene file '%s' generated successfully.\n", pov_filename);

    // --- Call POV-Ray to Render ---
    char command[512];
    memset(command,0,512);
    // Adjust the path to povray.exe if it's not in your system's PATH
    // On Linux/macOS: "povray -W1200 -H400 +A +R -O%s %s"
    // On Windows: "pvengine /RENDER /EXIT /W1200 /H400 /O%s %s"
    // I'll use a generic command that works on many systems if 'povray' is in PATH.
    // Make sure to replace `povray` with `pvengine` if you're on Windows and it's not in your PATH.
    sprintf(command, "/usr/bin/povray -W1800 -H300 -O%s %s", output_image_filename, pov_filename);

    printf("Executing POV-Ray command: %s\n", command);
    int result = system(command);

    if (result == 0) {
        printf("POV-Ray rendering successful. Image saved to '%s'\n", output_image_filename);
    } else {
        fprintf(stderr, "Error: POV-Ray rendering failed with code %d. Make sure POV-Ray is installed and in your system's PATH.\n", result);
        fprintf(stderr, "You might need to adjust the 'povray' command to 'pvengine' or specify its full path.\n");
    }

 return(0);
}







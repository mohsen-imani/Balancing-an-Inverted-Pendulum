
/*****************************************************************************/
/* File:        Reinf_learn.c                                                */
/* Description: Learning for Cart Pole System                                */
/* Author:                                                                   */
/* Date:                                                                     */
/* Modifications :                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Pole_sim.h"
#include "misc.h" 
#include "Reinf_learn.h"
#include <stdlib.h>

//***********************************MOHSEN***********************************
#define STATES  81
#define ACTIONS  2
#define MAX_TRAIL 1000
#define EPSILON .001

#define ALPHA 0.1
#define GAMMA 0.95
#define LAMDA 0.1
long double Q[STATES][ACTIONS];
long double e[STATES][ACTIONS];

int s_prev; // previous state
int s_current; // current state
int a_current = 0; // current action
int a_prev = 0; // previous action
int flag = 0; // a flag for determining when we are in the reset case

long double maxQ(int s_next);
int reward (int failure);
double random_number();
int  get_action(int s, double f[3], int invert, int *explore);
int  action_star(int s);
//***********************************MOHSEN***********************************
// This function initialize the Q  and eligibility matrix 
void initialize(){
int i ;
int j ;
    for(i = 0; i <= (STATES - 1); i++){
        for(j = 0; j <= (ACTIONS - 1); j++){
            Q[i][j] = 0.0;
            e[i][j] = 0.0;
		} // j
	} // i
}
/*****************************************************************************/
/* Multi-dimension decoder according to the discretization in Anderson, 1989 */
/* Input argument indicates the dimension to discretize (0 = X, 1 = Y, 2 = Z)*/
/* and the data structure for the pole-cart system.                          */
/* Return value is the state number (a total of 81 states per dimension).    */
/* Also computes the invert flag used to reduce states space using symmetry, */
/* and what is considered failure and writes it into  fail                   */
/* Input Variables:                                                          */
/*                  axis : the dimension to be encoded                       */
/*                  pole : the data structure of the cart-pole system        */
/* Output Variables:                                                         */
/*                  invert : inversion flag indicating the use of symmetry   */
/*                  fail   : fail flag indicates what is considered failure  */
/* Return Value: Numeric discretized state index for this dimension          */
/*****************************************************************************/


int Decoder3Dn(axis, pole, invert, fail)
int axis;
Polesys *pole;
int *invert;
int *fail;
{
  int pos, vel, ang, avel;
  static double pos_val[4] = {-1.5, -0.5, 0.5, 1.5};
  static double vel_val[4] = {-9999999999.9, -0.5, 0.5, 9999999999.9};
  static double ang_val[7] = {-0.20943951, -0.10471976, -0.017453293, 0.0, 
				        0.017453293, 0.10471976, 0.20943951};
  static double avel_val[4] = {-9999999999.9, -0.87266463, 0.87266463, 
				 9999999999.9};
	
  pos = -1;
  while ((pos < 3) && (pos_val[pos+1] < pole->pos[axis])) 
    ++pos;
  vel = -1;
  while ((vel < 3) && (vel_val[vel+1] < pole->vel[axis])) 
    ++vel;
  if (axis < 2) {
    ang = -1;
    while ((ang < 6) && (ang_val[ang+1] < (pole->theta[1-axis]
					   -(double)(axis)*0.5*M_PI))) 
      ++ang;
    avel = -1;
    while ((avel < 3) && (avel_val[avel+1] < pole->theta_dot[1-axis])) 
      ++avel;
  }
  else {
    ang = -1;
    while ((ang < 6) && (ang_val[ang+1] < MAX(fabs(pole->theta[1]), 
					      fabs(pole->theta[0]-0.5*M_PI)))) 
      ++ang;
    avel = -1;
    while ((avel < 3) && (avel_val[avel+1] < 
			  MAX(SIGN(pole->theta[1])*pole->theta_dot[1],
			      SIGN(pole->theta[0])*pole->theta_dot[0]))) 
      ++avel;
  }
    
  // Sets fail, i.e. if the trial should be considered to have ended based on 
  // this dimension
  *fail = ((pos == -1) || (pos == 3) || (vel == -1) || (vel == 3) || (ang == -1)
 	  || (ang == 6) || (avel == -1) || (avel == 3));
   
  // Use symmetry to reduce the number of states

  if (!(*fail))
    {
      *invert = 0;
      if (ang > 2)
	     {
		   *invert = 1;
		   ang = 5-ang;
		   pos = 2-pos;
		   vel = 2-vel;
		   avel = 2-avel;
		 }
	  return(pos + 3*vel + 9*ang + 27*avel);
    }
	// Failed situations are not part of the state space
	return(-1);
}
       
       


/*****************************************************************************/
/* Main learning function. Takes the information of the system from   pole   */
/* and the   reset   flag which indicates that this is the first state in a  */
/* new  trial. The action to take in the next time step is written into the  */
/* force    vector and then applied by the simulator to the 3D cart. Also    */
/* returned is the information whether the trial should be ended using the   */
/* fail  flag  and a counter of the number of exploration actions that have  */
/* been taken within this trial by incrementing the  explore  counter every  */
/* an exploration action is taken.                                           */
/* Input Variables:                                                          */
/*                  pole  : the data structure of the cart-pole system       */
/*                  reset : flag indicating that this is the first step in a */
/*                          new trial (i.e. unrelated to the previous state  */
/*                  explore : the number of exploration akitions taken in    */
/*                          this trial prior to this time step               */
/* Output Variables:                                                         */
/*                  force : force vector to be applied to the cart in the    */
/*                          next time step (corresponding to the action taken*/
/*                  fail  : flag indicating whether a new trial should be    */
/*                          started in the next time step                    */
/*                  explore : the number of exploration taken in this trial  */
/*                            including this time step (increase by one if   */
/*                            exploration action was taken)                  */
/*****************************************************************************/




void pole_learn(pole, reset, force, fail, explore)
Polesys *pole;
int reset;
double force[3];
int *fail;
int *explore;
{

  int invert;
  int a_star;
  int i, j; 
  long double delta;
// set the previous state and force and actions in the case of Reset
if (reset == 1){
   /*init_pole();

   s_prev = Decoder3Dn(0, pole, &invert, fail);
   force[1] = 0.0;
   force[2] = 0.0;     
   force[0] = 0.0;
   a_prev = ((force[0] - F_X) > -0.001 && (force[0] - F_X) < 0.001)? 1 : 0;
   s_current = s_prev ;
   a_prev = a_current = 0;*/
   flag = 1;
   return;
} 
// updating Q and eligibility
  if (flag == 1) { // Here we are in the initial state
     s_prev = Decoder3Dn(0, pole, &invert, fail);
     a_prev = get_action(s_prev, force, invert, explore);
     flag = 0;
     return;
   }
  if (flag == 0){ // Here we have passed the initial state and we can update the Q
     s_current = Decoder3Dn(0, pole, &invert, fail);
     s_current = (*fail == 0)? s_current: s_prev;
     a_current = (*fail == 0)? get_action(s_current, force, invert, explore): a_prev;
     a_star = action_star(s_current);

     delta = ((long double)reward(*fail) + GAMMA * Q[s_current][a_star] - Q[s_prev][a_prev]);
     e[s_prev][a_prev] += 1;

     for (i = 0; i < STATES; i++){
         for (j = 0; j < ACTIONS; j++) {
            Q[i][j] += delta * ALPHA * e[i][j];
            e[i][j] = (a_current == a_star) ? (GAMMA * LAMDA * e[i][j]) : 0.0;

         }
     }
     a_prev = a_current;
     s_prev = s_current;
   }
  return; 
}	
 


long double maxQ(int s_next){
    return (Q[s_next][0] > Q[s_next][1]? Q[s_next][0] : Q[s_next][1]);
    //printf("maxQ[%d] = %Lf\n", s_next, out);    
} 

int reward (int failure){
    return (failure == 1 ? -1 : 0);
}

double random_number()
{
    return ((double)rand() / (double)RAND_MAX);;
}

int  get_action(int s, double f[3], int invert, int *explore){
    long double q_L = Q[s][0];  
    long double q_R = Q[s][1];
    f[1] = 0.0;
    f[2] = 0.0;
    if (random_number() > EPSILON) {
       f[0] = (invert == 0)?((q_L > q_R)? -(F_X):F_X):-((q_L > q_R)? -(F_X):F_X);
       return ((q_L > q_R)? 0 : 1); 
    } else {
       int r = (rand() % 100+1);
       f[0] = ( r > 50? F_X:-(F_X));
       (*explore)++;

       return (r > 50? 1 : 0);
     }
}

// this function returns the best action (a*)
int  action_star(int s){
     return ((Q[s][0] > Q[s][1])? 0 : 1);
}

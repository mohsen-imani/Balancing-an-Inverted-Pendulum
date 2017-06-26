
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
#define STATES  162
#define ACTIONS  2
#define MAX_TRAIL 1000
#define EPSILON .9999
#define ALPHA 0.1
#define GAMMA 0.95
#define BETA 0.2
#define TEMP 0.001
#define LAMDA 0.5

double value[STATES];
double pref[STATES][ACTIONS];
long double e_actor[STATES][ACTIONS];
long double Q[STATES][ACTIONS];
long double e_critic[STATES];


long double maxQ(int s_next, int f);
int reward (int failure);
double random_number();
int  get_action(int s, double f[3], int *invert, int *explore);
int  get_action_AC(int s, double f[3], int *invert, int *explore);
void updateV_value(int s);
double pr(int s, int a);

void initialize(){
int i ;
int j ;
    for(i = 0; i <= (STATES - 1); i++){
        for(j = 0; j <= (ACTIONS - 1); j++){
            pref[i][j] = 0.0;
            e_actor[i][j] = 0.0;            
		} // j
        value[i] = 0.0;
        e_critic[i] = 0.0;
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
	  return(pos + 3*vel + 9*ang + 27*avel);
      //if (ang > 2)
	    /* {
		   *invert = 1;
		   ang = 5-ang;
		   pos = 2-pos;
		   vel = 2-vel;
		   avel = 2-avel;
		 }*/

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

//*********************MOHSEN***************************
int s_prev; // previous state
int s_current; // current state
int a_current = 0; // current action
int a_prev = 0; // previous action
int flag = 0;

//*********************MOHSEN***************************


void pole_learn(pole, reset, force, fail, explore)
Polesys *pole;
int reset;
double force[3];
int *fail;
int *explore;
{

  int invert,i ,j; 
  long double delta; 

if (reset == 1){
   init_pole();
   flag = 1;
   return;
} 
// updating Q and eligibility
  if (flag == 1) { // Here we are in the initial state
     s_prev = Decoder3Dn(0, pole, &invert, fail);
     a_prev = get_action_AC(s_prev, force, invert, explore);
     flag = 0;
     return;
   }
  //printf("/------------------------------------------------/");
  s_current = Decoder3Dn(0, pole, &invert, fail);
  a_current = (*fail == 0)? get_action_AC(s_current, force, invert, explore): a_prev;
  s_current = (*fail == 0)? s_current:s_prev;
  
if (flag == 0){ 
  
  delta =  (long double)reward(*fail) + GAMMA *  value[s_current] - value[s_prev];
  e_critic[s_prev] += 1;


  pref[s_prev][a_prev] =  pref[s_prev][a_prev] + ALPHA * delta * e_actor[s_prev][a_prev];
  e_actor[s_prev][a_prev] += 1;
  
  for (i = 0; i < STATES; i++){
      value[i] += ALPHA * delta * e_critic[i];
      e_critic[i] = GAMMA * LAMDA * e_critic[i];
      for (j = 0; j < ACTIONS; j++) {
          //Q[i][j] += delta * ALPHA * e[i][j];
          e_actor[i][j] = GAMMA * LAMDA * e_actor[i][j];
      }
  } 
  }

  if (s_current != -1){ 
     a_prev = a_current;
     s_prev = s_current;
     flag = 0;
   }
  return; 
}	
 


long double maxQ(int s_next, int f){
long double out;
if (f ==0){
    out = (Q[s_next][0] > Q[s_next][1]? Q[s_next][0] : Q[s_next][1]);
    printf("maxQ[%d] = %Lf\n", s_next, out);    
    return out;
} else {
    printf("maxQ[%d] = 0\n", s_next);
    return 0.0;
}
} 

int reward (int failure){
    return (failure == 1 ? -1 : 0);
}

double random_number()
{
    return ((double)(rand() % 10000+1) / 10000.0);
}

int  get_action(int s, double f[3], int *invert, int *explore){
    long double q_L = Q[s][0];  
    long double q_R = Q[s][1];
    //printf("get_action s = %d q_L = %Lf q_R = %Lf",s,q_L,q_R);
    f[1] = 0.0;
    f[2] = 0.0;
    if (random_number() < EPSILON) {
       f[0] = ((q_L > q_R)? -(F_X):F_X);
    //printf("  force %f action %d\n",f[0],((q_L > q_R)? 0 : 1));
       return ((q_L > q_R)? 0 : 1); 
    } else {
       int r = (rand() % 100+1);
       f[0] = ( r > 50? F_X:-(F_X));
       (*explore)++;
    //printf("  force %f action %d\n",f[0],(r > 50? 1 : 0));
       return (r > 50? 1 : 0);
     }
}

int  get_action_AC(int s, double f[3], int *invert, int *explore){ 
    double rndm = random_number();
    //printf("get_action s = %d q_L = %Lf q_R = %Lf",s,q_L,q_R);
    f[1] = 0.0;
    f[2] = 0.0;
       f[0] = ((rndm < pr(s,0))? -(F_X):F_X);
    //printf("  force %f action %d\n",f[0],((q_L > q_R)? 0 : 1));
       return ((rndm < pr(s,0))? 0 : 1); 
   /* } else {
       int r = (rand() % 100+1);
       f[0] = ( r > 50? F_X:-(F_X));
       (*explore)++;
    //printf("  force %f action %d\n",f[0],(r > 50? 1 : 0));
       return (r > 50? 1 : 0);
     }*/
}

void updateV_value(int s){
    value[s] = pr(s,0)*Q[s][0] +pr(s,0)*Q[s][1];
    return;
    }
// This function returns the probability of selection action a in state s 
double pr(int s, int a){
   return (pow (2.71828, pref[s][a]/TEMP)/(pow (2.71828, pref[s][0]/TEMP) + pow (2.71828, pref[s][1]/TEMP)));
} 
      
